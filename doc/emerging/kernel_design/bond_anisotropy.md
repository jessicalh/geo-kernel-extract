Historical -- not current truth; see doc/emerging/kernel_design/mcconnell_spec.md and doc/emerging/CONTROLLING_SPEC.md.

# Bond magnetic anisotropy — kernel design

What this is: a design note for **one input feature** to an e3nn equivariant
shielding predictor — the through-space shielding a nucleus feels from the
anisotropic magnetic susceptibility of an individual *localized bond* nearby
(C=O, C≡N, C–C, C=C, etc.), the "bond anisotropy" or "neighbour anisotropy"
contribution evaluated bond-by-bond in the ApSimon / McConnell point-dipole
tradition. The brief is "what *should* we compute, the defensible standard way",
not a verdict on our code. Web-grounded first, code-read second.

**Read `mcconnell.md` first.** This note is mostly a pointer into it. Bond
anisotropy and the McConnell neighbour-anisotropy feature are **the same
machine** — the same Δχ-tensor ⊗ dipole-propagator → `2e` construction — differing
only in *who the source is* (a single chemical bond, with a per-bond Δχ) and in
a few bond-specific conventions (where on the bond the dipole sits; axial vs
rhombic per bond type; which bonds are already covered by ring/McConnell so we
don't double-count). The bulk of the physics, the irrep choice, and the honest
ceilings are stated once in `mcconnell.md` and not repeated here.

Terms are defined as they appear; "Δχ", "source", "target", "T0/T1/T2",
"`0e`/`1o`/`2e`" carry the same meanings as in `mcconnell.md`.

---

## 1. The effect, and how the literature says to compute it

### 1.1 The physical picture

Every bond is a little anisotropic magnetic susceptibility. In the field B₀ its
electrons circulate, and because the bond responds more strongly along some
local axes than others (e.g. more along a C≡C axis than across it), the induced
moment is misaligned with B₀. That misaligned moment is a magnetic dipole; its
stray field at a nearby nucleus shifts the resonance. The "anisotropy of the
neighbouring bond" is exactly the through-space deshielding/shielding cone you
see around C=O, C=C, C≡C, and (weakly) C–C bonds, and the reason aldehyde
protons sit at 9–10 ppm, vinyl protons at 4.5–6.5 ppm, and alkyne protons
anomalously upfield at 2–3 ppm
([Chem LibreTexts 14.8, diamagnetic anisotropy](https://chem.libretexts.org/Bookshelves/Organic_Chemistry/Map%3A_Organic_Chemistry_(Bruice)/14%3A_NMR_Spectroscopy/14.08%3A_Diamagnetic_Anisotropy)).
Physically it is the *same* effect as the aromatic ring current and the McConnell
group-anisotropy term — a localized current loop's dipole field — just applied to
a two-centre bond rather than a delocalised ring or a whole functional group
([RSC Adv. 2014, C3RA45512B](https://pubs.rsc.org/en/content/articlehtml/2014/ra/c3ra45512b)).

### 1.2 The standard computation: McConnell point dipole, per bond

The field-standard recipe is McConnell's 1957 point-dipole model applied to a
single bond ([McConnell, *J. Chem. Phys.* 27, 226](https://doi.org/10.1063/1.1743676)).
Collapse the bond's induced currents to a point dipole sitting on the bond axis,
evaluate its dipole field at the target. For an **axially symmetric** bond (one
Δχ describing it) the orientation-averaged scalar shift is the familiar form
(here in the protein-predictor notation,
[Abraham/Modgraph SCS7, *C–C bond anisotropy*](https://www.modgraph.co.uk/Downloads/SCS7.pdf)):

```
δ  =  Δχ_bond · (1 − 3 cos²φ) / (3 · L₀ · R³)
```

where R is the dipole→nucleus distance, φ the angle from the bond axis, Δχ_bond
the molar anisotropy of that bond category, and L₀ Avogadro's number (the molar
prefactor). This is the **same (1−3cos²θ)/r³ object as McConnell and the
point-dipole ring** — its full, un-averaged form is the rank-2 dipole propagator
**D(r) = (3 r̂⊗r̂ − I)/r³** contracted with the bond's Δχ tensor, and the scalar
above is its axial, orientation-averaged projection. See `mcconnell.md` §1.2–1.4
for the full-tensor and spherical-harmonic (Σ χ_m Y₂^m) statements; they apply
verbatim, the only change being that the susceptibility tensor is a single bond's
rather than a group's.

### 1.3 The two genuine bond-specific wrinkles

Two things are specific to the **per-bond** treatment and are not fully spelled
out in `mcconnell.md`:

1. **Axial (one Δχ) vs non-axial (two Δχ) by bond type.** The literature is
   explicit and consistent: bonds with cylindrical symmetry — the triple bond
   C≡C / C≡N — are described by **one** anisotropy Δχ along the bond axis; but
   non-cylindrical groups, paradigmatically the **carbonyl C=O**, need **two
   distinct anisotropies**, Δχ_parallel and Δχ_perpendicular (i.e. χ has three
   inequivalent principal values, an axial *and* a rhombic part), because the
   carbonyl is flat, not rod-symmetric
   ([Abraham Part 13, *ketones — anisotropy and electric field of C=O*](https://pubs.rsc.org/en/content/articlehtml/1999/p2/a808908f);
   [RSC Adv. 2014, C3RA45512B](https://pubs.rsc.org/en/content/articlehtml/2014/ra/c3ra45512b)).
   So a defensible bond-anisotropy feature is **axial for C≡C/C≡N/C–C and rhombic
   for C=O/C=C** — the rhombicity is *tabulated*, not invented, for exactly the
   bonds where it matters. This is the one place a per-bond feature carries more
   structure than the axial-only default in `mcconnell.md` §2.1.

2. **Where the dipole sits on the bond.** McConnell collapses the bond to a
   *point*; the point's position along the bond axis is a model choice that
   matters at short range. ApSimon, and Abraham after him, do not put it at one
   fixed end — for the carbonyl Abraham fits the dipole origin ~1.0 Å from the
   carbon along the C=O axis, with separate effective parallel/perpendicular
   components placed there
   ([Abraham/Modgraph, *aliphatic ketones*](http://www.modgraph.co.uk/Downloads/Aliphatic%20ketones.pdf)).
   Our existing McConnell code uses the **bond midpoint** as the origin (§4); the
   defensible options are midpoint (geometric, parameter-free) or a fitted
   on-axis offset (matches Abraham/ApSimon but adds a per-category parameter). For
   an equivariant feature where the magnitude is going to be calibrated anyway
   (§2), the midpoint is the cleaner default and the offset is a refinement, not a
   necessity.

### 1.4 ApSimon: the same model, a sharper source treatment

ApSimon's contribution (the name in the brief) is **not a different physics** — it
is the McConnell point-dipole model done more carefully for bonds: a
**double-summation** evaluation of the susceptibility over the bond's two centres
and a fit of *two* anisotropy components for non-axial bonds rather than one. Its
headline finding for the C=C double bond is that the naive single-cone picture
"requires substantial modification" — the deshielding is confined to the bond
*ends*, with shielding elsewhere in and above the plane — precisely the rhombic
(two-Δχ) correction of point (1) above
([RSC Adv. 2014, C3RA45512B](https://pubs.rsc.org/en/content/articlehtml/2014/ra/c3ra45512b);
locally held: Martin 2000, *graphical model for proton deshielding over C=C*,
`references/` — the modern restatement of the ApSimon alkene result). For our
purposes ApSimon ⇒ "carry the rhombic Δχ for C=O and C=C", nothing more exotic.

### 1.5 Where the Δχ values come from (and the double-count trap)

Per-bond Δχ is tabulated, not derived. The canonical compilations and re-fits are
**Zürcher** (steroid methyl probes, the pioneering set), **ApSimon**,
**Schneider**, **Williamson–Asakura** (the protein peptide-bond lineage: paired
C=O and C–N constants in the amide plane), and **Abraham**'s repeated re-fits
(ketones, aliphatic/aromatic amides). The project already holds these as a
research note — **`MCCONNELL_DCHI_LITERATURE.md`** — with converted molar-cgs
values, the molar→ppm prefactor, paired C=O/C–N peptide constants, and a per-row
confidence column. That note is the bond-anisotropy Δχ source; this design does
not duplicate it. Two facts from it that are load-bearing here:

- The values **disagree across sources by ~2–5×** for the same bond (carbonyl
  C=O spans roughly +2 to +24 ×10⁻⁶ cm³ mol⁻¹ across Williamson / Abraham /
  ApSimon / Schneider / Zürcher). The scatter *is* the honest statement of how
  well the parameter is known; it argues for Δχ as a **calibrated coefficient**,
  not a constant.
- **Aromatic bonds must be zeroed if a ring-current kernel is active.** The large
  aromatic susceptibility anisotropy *is* the π ring current re-expressed as a
  bond property; counting both double-counts the same physics. `MCCONNELL_DCHI_LITERATURE.md`
  flags this explicitly and the project's ring kernels already carry the
  per-ring intensities.

---

## 2. The defensible, sane recommendation

**Bond anisotropy collapses into the McConnell feature.** Compute it as the
per-bond instance of the McConnell construction in `mcconnell.md` §2: for each
bond within a KD-tree cutoff, build the bond's susceptibility tensor and contract
it with the dipole propagator **D(r) = (3 r̂⊗r̂ − I)/r³**, accumulate per bond
category, and emit the irreducible pieces. Do **not** stand up a parallel
calculator and do **not** start from (1−3cos²φ)/r³. Concretely, departing from
`mcconnell.md` only where the per-bond physics genuinely differs:

1. **Source tensor — axial or rhombic *by bond category*.** Build
   `Δχ_tensor(bond) = Δχ_ax · (ŝ ⊗ ŝ − I/3)  +  Δχ_rh · (rhombic part)` where ŝ
   is the bond axis. Use **axial only** (Δχ_rh = 0) for C≡C, C≡N, C–C, S–S; carry
   the **tabulated rhombic** Δχ_rh for **C=O and C=C**, where the literature
   actually supplies two components (§1.3). Where no rhombic value is tabulated,
   Δχ_rh = 0 and say so — that is honest absence, not an approximation chosen for
   convenience.

2. **Geometry emitted unscaled (Å⁻³); Δχ is a calibrated per-category
   coefficient.** Exactly as McConnell, ring, and the project's
   literature-scaling pattern. The ~2–5× source scatter (§1.5) becomes a fitted
   coefficient with a reported band, not a hard-coded guess. Pull candidate
   values from `MCCONNELL_DCHI_LITERATURE.md`.

3. **Irreps — primary `2e`, plus the parity caveat.** The clean, single-ℓ form is
   the symmetric-traceless propagator-on-Δχ tensor decomposed to its five real
   spherical components, emitted as **`2e`**, per bond category. The orientation-
   averaged scalar (1−3cos²φ)/r³ is recoverable as the `0e` projection and is a
   *sanity baseline*, not the primary feature for an equivariant model. Do not
   synthesise a `1o` from cross terms unless a named asymmetry justifies it. All
   of §2.2 and §3 of `mcconnell.md` apply unchanged.

4. **Per-category channels, with double-count discipline.** Emit per bond
   category (C=O, C–N, C≡C/C≡N if present, C–C, S–S) so the model carries a
   different Δχ per chemistry. **Zero aromatic bonds while the ring kernel is
   active** (§1.5). And be explicit about the overlap with the existing McConnell
   group term: McConnell already iterates these same bonds (§4), so a separate
   "bond anisotropy" emit would re-sum bonds McConnell has summed — the right move
   is to **extend the single McConnell calculator** (richer per-bond Δχ, rhombic
   C=O/C=C, the S–S and C–C categories it does not yet weight), not to add a
   second kernel over the same bonds.

**Dipole origin (§1.3 point 2):** default to the bond **midpoint** (parameter-
free, what our code does); treat the Abraham/ApSimon fitted on-axis offset as an
optional refinement, not a requirement, since the magnitude is calibrated
downstream anyway.

---

## 3. Where the physics is clean and where it is not

**Clean:**
- The **rank-2 (`2e`) angular structure** is exact and parameter-free given
  geometry — the dipole propagator is rigorous classical EM, same as McConnell
  and the point-dipole ring.
- The **axial-vs-rhombic distinction by bond type** is itself well-established and
  tabulated (one Δχ for C≡C/C≡N, two for C=O/C=C). Carrying it is standard, not
  novel.

**Not clean (footnotes, not bugs):**
- **Δχ magnitude.** Tabulated values disagree ~2–5× across Zürcher / ApSimon /
  Schneider / Williamson / Abraham (§1.5). Response: calibrated per-category
  coefficient with a band, not a constant.
- **Rhombicity coverage.** Tabulated for C=O and C=C; absent for most single
  bonds. Δχ_rh = 0 there is honest absence.
- **Point vs distributed / dipole-origin placement.** The point-dipole collapse
  is a far-field model; for a two-centre bond at the 2.5–4 Å contacts that matter
  in proteins it is in its weakest regime, and *where* on the bond axis the point
  sits (midpoint vs a fitted offset) measurably moves the near field. We disclose
  this rather than engineer it away; a distributed/GIAO-TSNMR treatment is a
  different, heavier method, off the table for a standard non-novel feature
  (`mcconnell.md` §1.5).
- **Overlap / double-counting.** This is the *defining* not-clean issue for this
  kernel. Bond anisotropy is the same machine as McConnell over the same bonds,
  and aromatic bond anisotropy is the same physics as the ring current. Both
  overlaps must be handled by *extending the one McConnell calculator and zeroing
  aromatic bonds*, or the model sees the same field two or three times.
- **Parity of emitted irreps.** A susceptibility tensor is even; an outer-product-
  of-displacements tensor is not. Declare the parity correctly to e3nn
  (`mcconnell.md` §3).

**One-line summary:** bond anisotropy is the McConnell neighbour-anisotropy
feature with a per-bond Δχ table (and a tabulated rhombic part for C=O/C=C); the
angular form is the exact `2e` dipole propagator, the magnitude and the point-
source/origin assumptions are the soft parts, and the central design obligation
is *not to double-count* the bonds McConnell and the ring kernel already cover.

---

## 4. Neutral note: current compute vs this recommendation

Descriptive only, no verdict. There is **no separate "bond anisotropy"
calculator** — the library's `src/McConnellResult.{h,cpp}` *is* the per-bond
anisotropy kernel. Its header states it plainly: "For each bond within the
configured cutoff (default 10 Å), computes the full McConnell tensor kernel."
It loops over bonds by `BondCategory` (`Types.h`: PeptideCO, PeptideCN,
BackboneOther, SidechainCO, Aromatic, Disulfide, SidechainOther, Unknown), takes
the **bond midpoint** as the dipole origin and the **bond direction** as the
axis, and forms the three-term tensor `M = 9 cosθ (d̂⊗b̂) − 3 (b̂⊗b̂) − (3 d̂⊗d̂ − I)`
stored as `M/r³` (Å⁻³), keeping the symmetric-traceless dipolar kernel K
separately and the scalar f = (3cos²θ−1)/r³ as K's double-contraction with the
bond axis. It accumulates per-category sums (CO/CN/sidechain/aromatic) and emits
the 9-component decomposition, per-category symmetrised-traceless `2e` blocks,
and the category scalar sums. Susceptibility scaling is **not** applied by the
producer; the geometry is emitted unscaled, and the per-bond Δχ candidate values
live in the rediscover-side note `MCCONNELL_DCHI_LITERATURE.md` (Williamson–
Asakura, Abraham, ApSimon, Schneider, Zürcher, with the molar→ppm prefactor and
the aromatic double-count flag). Relative to §2, the descriptive differences are:
(a) the current source model is **axial-by-construction** (one Δχ per category,
applied via the geometry; no tabulated **rhombic** Δχ_rh carried for C=O/C=C —
the wrinkle of §1.3); (b) the **Disulfide** and single C–C categories exist in
the enum but are not weighted as distinct anisotropy sources; (c) the primary
emitted tensor is the full three-term M (carrying T0+T1+T2 together) alongside
the separate clean K, rather than the clean `2e` propagator-on-Δχ being named the
single production feature (the same observation `mcconnell.md` §4 makes). The
dipole origin is the midpoint, which matches §2's recommended default. These are
descriptions, not conclusions.

---

## Sources

- McConnell, *J. Chem. Phys.* 27, 226 (1957) — point-dipole model.
  [DOI 10.1063/1.1743676](https://doi.org/10.1063/1.1743676)
- RSC Adv. 2014, C3RA45512B — neighbour/bond anisotropy review; axial vs
  non-axial (carbonyl two-component) bonds; ApSimon C=C re-analysis; magic angle.
  [C3RA45512B](https://pubs.rsc.org/en/content/articlehtml/2014/ra/c3ra45512b)
- Abraham et al., *Proton chemical shifts Part 13 — ketones, C=O anisotropy and
  electric field* — two-component carbonyl Δχ, dipole origin on the C=O axis.
  [Perkin Trans. 2 1999, A808908F](https://pubs.rsc.org/en/content/articlehtml/1999/p2/a808908f),
  [aliphatic ketones PDF](http://www.modgraph.co.uk/Downloads/Aliphatic%20ketones.pdf)
- Abraham/Modgraph SCS7 — C–C and C–H single-bond anisotropy, the
  δ = Δχ(1−3cos²φ)/3L₀R³ bond form. [SCS7 PDF](https://www.modgraph.co.uk/Downloads/SCS7.pdf)
- Chem LibreTexts 14.8 — diamagnetic anisotropy of C=C, C≡C, C=O (the
  observable shift signatures). [LibreTexts 14.8](https://chem.libretexts.org/Bookshelves/Organic_Chemistry/Map%3A_Organic_Chemistry_(Bruice)/14%3A_NMR_Spectroscopy/14.08%3A_Diamagnetic_Anisotropy)
- Martin 2000, *graphical model for proton deshielding over a C=C double bond* —
  modern restatement of the ApSimon alkene result (held locally in `references/`).
  [MDPI IJMS 1(4):84](https://mdpi.com/1422-0067/1/4/84/htm)
- Case, *Chemical shifts in biomolecules* — full-tensor neighbour form, isotropic-χ
  cancellation. [PMC3877577](https://pmc.ncbi.nlm.nih.gov/articles/PMC3877577/)
- e3nn docs; Tensor Field Networks (arXiv:1802.08219) — rank-2 → `0e ⊕ 1 ⊕ 2e`.
  [e3nn](https://docs.e3nn.org/), [arXiv:1802.08219](https://arxiv.org/abs/1802.08219)
- Project-internal (not literature): `mcconnell.md` (the parent design),
  `MCCONNELL_DCHI_LITERATURE.md` (per-bond Δχ tables + prefactor + double-count
  flag), `src/McConnellResult.{h,cpp}`, `src/Types.h` (`BondCategory`).
