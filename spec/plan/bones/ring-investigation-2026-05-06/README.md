# Ring-investigation framing memo (2026-05-06)

**Purpose.** Before drafting Bundle C (substrate-driven ring construction
to replace `Protein::DetectAromaticRings`, plus Pro pyrrolidine
adoption), we need empirical knowledge of every calculator's actual
ring interaction. Bundle C cannot be planned on speculation
("calculators may filter," "some emit choices") — those phrases were
the antipattern that triggered this investigation.

This memo is the framing document every investigation agent reads
first. It establishes the substrate vocabulary, the typed-identity
discipline, the science background each agent needs, the per-report
template, and the investigation process.

## Status anchor

Repository: `/shared/2026Thesis/nmr-shielding/`. As of 2026-05-06:

- **Bundle B substrate runtime population is landed.**
  `LegacyAmberTopology` carries per-atom `AtomSemanticTable` populated
  via typed `LookupBy / LookupCap / ApplyCapDelta` at
  `Protein::FinalizeConstruction`. Substrate composition reads
  `Atom.pdb_atom_name` once via `ParseAtomName` at the boundary;
  thereafter typed identity flows.
- **Naming-application architecture is landed.** `NamingRegistry.{h,cpp}`
  holds `NamingApplicator` with rule-application architecture:
  `NamingSource`-tagged rules each with `Applies()` predicate +
  `Output()` action; per-atom transient application map collected via
  `Collect()`; explicit `Resolve()` method body documents per-case
  decisions; post-resolution validator enforces canonical-output
  invariant; deletion-variant deny overlay in `IsCanonical`. Architecture
  is `spec/plan/naming-applicator-architecture-sketch-2026-05-06.md`.
- **Pre-Bundle-C state.** `Protein::DetectAromaticRings` is still
  string-matched against `AminoAcidType::rings[].atom_names` const-
  char-pointer arrays. This is what Bundle C replaces. The substrate
  has all the typed information needed for substrate-driven ring
  construction; Bundle C wires it.
- **Pro pyrrolidine adoption** is part of Bundle C per the user's
  "no blur" rule on substrate-encoded chemistry. `RingTypeIndex`
  gains `ProPyrrolidine`; new `ProPyrrolidineRing` class; Pro
  residues materialise a `Ring` object with `Intensity = 0`,
  `Aromaticity = None`, `RingSize = 5`, `NitrogenCount = 1`.
- **`AminoAcidType::rings[].atom_names`** removal is part of Bundle C
  (after investigation confirms only `DetectAromaticRings` reads it).
  `ChiAngleDef::atoms[]` (parallel string surface for chi angles)
  STAYS — that's audit Hotspot 2, separate slice.

## Why investigation precedes Bundle C

The Bundle C plan as drafted before this investigation made
unverified claims of the form "calculators may filter Pro rings,"
"some may iterate inside KernelFilterSet," "some may have hardcoded
[8] arrays." These speculations were architectural fragility:
Bundle C cannot assume calculator behaviour without inspection.

In particular, the author misunderstood `GeometryChoice` as side-
effect noise; it is the project's deliberate audit surface. Per
PATTERNS.md "When implementing a calculator":

> You must record GeometryChoices via GeometryChoiceBuilder during
> Compute(). Every inclusion, exclusion, and triggered event gets a
> Record() call inside a lambda. Attach the entities involved
> (atoms, rings, bonds) with their roles and outcomes. Add named
> numbers (distance, intensity, threshold) with units.

A calculator iterating a Pro ring and skipping it because of
`Aromaticity = None` ENRICHES the choice record with documentary
value (future calculators can ask "what was considered? what was
rejected? why?"). Pro rings flowing through current calculators is
an information gain, not a defect — provided each calculator's
filter logic correctly produces zero contribution under
`Intensity = 0`.

The investigation establishes empirically what each calculator
does.

## Science grounding — why ring chemistry is load-bearing for this thesis

This section is required reading. Each agent uses it to ground their
calculator's behaviour in what the thesis actually claims, so the
investigation report can answer "what is this calculator's
contribution to the thesis result?" instead of just "what does this
calculator compute?". No handwaving — every claim below is concrete
with citations or pointers to in-repo verifiable values.

### 1. The thesis's primary result is the T2 angular residual

This is not a slogan. From `spec/MATHS_GOALS.md` and `PATTERNS.md`
"T2 Completeness — This Is Not Optional":

A geometric kernel evaluated at an atom produces a 3×3 tensor.
Decomposed into spherical harmonics: T0 (isotropic, 1 component),
T1 (antisymmetric pseudovector, 3 components), T2 (symmetric
traceless, 5 components). T0 is the classical chemical shift —
every NMR predictor since Pople (1956) computes some version. T2 is
what this thesis adds: the angular pattern that reveals where the
classical model breaks down.

The thesis calibrates classical-kernel parameters against DFT
WT-ALA delta tensors and reports the T2 residual per atom-class —
the angular structure DFT computes that the classical sum cannot
match. For ring-current calculators specifically, the T2 shape
comes directly from ring vertex geometry through the kernel
integral. Ring atom ordering, ring normal direction, and ring
position labels are not display details — they ARE the angular
structure being measured.

Ring-current calculators that produce wrong T2 shape (because their
ring atom_indices are reverse-cyclic, because their normal sign is
flipped, or because they pool across ring positions) corrupt the
primary thesis output. Calibration appears to converge but the
fitted parameters absorb the error rather than reflect physics.

**Investigation implication.** Each per-calculator report's "Output
surface" section names every NPY column and ConformationAtom field
that carries T2 components. These are the surfaces the calibration
pipeline reads. Bundle C must preserve their bit-identity for non-Pro
rings; Pro rings (Intensity = 0) cannot perturb them.

### 2. Ring currents dominate aromatic-region NMR shielding

Concrete numbers from `PATTERNS.md` §22-§23 and the GEOMETRIC_KERNEL
_CATALOGUE:

- Phe ring with literature intensity I = -12.0 nA/T:
  - Probe atom 3 Å above ring face → σ = +1.40 ppm (shielded; correct
    sign verified analytically).
  - Probe atom 5 Å in-plane → σ = -0.16 ppm (deshielded; sign
    transitions correctly across the ring axial vs equatorial
    region).
- Backbone amide H typical chemical-shift range: 6-10 ppm.
- Ring-current shifts at H atoms 3-5 Å from aromatic ring centres:
  ±1-3 ppm (10-50% of the typical amide-H shift range for nearby
  protons).
- ¹³C ring-current shifts at backbone Cα: 0.5-2 ppm; at sidechain
  carbons within 5 Å of an aromatic ring: up to 5-8 ppm.

Reference: Case, D. A. (1995) "Calibration of ring-current models for
the calculation of protein chemical shifts," J. Biomol. NMR 6, 341-
346. This is the calibration paper from which the project's
literature ring-current intensities (`Ring::LiteratureIntensity()`
per RingTypeIndex) are taken.

For 25 typical proteins, ring-current contributions to ¹H shielding
near aromatic residues account for 20-40% of the total non-secondary-
structure shift signal (this fraction is what motivates the
calibration of geometric ring-current kernels against DFT;
non-classical effects + electronic polarization make up the residual
that the upstream e3nn model fits).

**Investigation implication.** A calculator that processes ring
chemistry produces a substantial fraction of the per-atom shielding
prediction. Errors in ring-atom identity, normal sign, or ring-type
parameterisation propagate to shielding predictions at every atom
within ~15 Å of an aromatic residue (the project's
`RING_CALC_CUTOFF_A = 15.0`). Calculator behaviour on Pro pyrrolidine
matters less for current calculators (Intensity = 0 → zero
contribution) but matters enormously for any future calculator that
would consume the Pro Ring (puckering, ring-flip, sterics).

### 3. Per-ring-type angular structure differs

From `Ring.h` subclass overrides plus literature:

```text
Ring type      I (nA/T)   ring_size   nitrogens   aromaticity
PheBenzene     -12.0      6           0           Full
TyrPhenol      -11.28     6           0           Full
TrpBenzene     -12.48     6           0           Full
TrpPyrrole     -6.72      5           1           Reduced
TrpPerimeter   -19.2      9           1           Full (synthesized)
HisImidazole   -5.16      5           2           Weak
HidImidazole   -5.16      5           2           Weak
HieImidazole   -5.16      5           2           Weak
ProPyrrolidine 0.0        5           1           None (saturated)
```

Per-ring-type T2 structure differs because:

- **Ring size determines vertex count** (6 vs 5 vs 9). Vertex count
  changes the angular sweep over which the line integral (Biot-
  Savart) or surface integral (Haigh-Mallion) accumulates.
- **Heteroatom count distorts the ring's electron density**.
  Imidazole (2 N) has different π density at carbons vs nitrogens
  than benzene (0 N). The point-dipole-at-centre model of
  RingSusceptibility approximates this via per-type Δχ; the
  vertex-walk Biot-Savart picks it up via per-vertex partial
  charges in some derivations.
- **Lobe offset (`JBLobeOffset`) differs** by ring type — 0.64 Å for
  6-rings, 0.50-0.52 Å for 5-rings, 0.60 Å for the TRP perimeter
  (intermediate). The two-loop Johnson-Bovey model places loops at
  ±lobe_offset from the ring plane; this offset is calibrated per
  ring type.

The T2 cosine similarity between calculators (PATTERNS.md §22):

```text
McConnell vs Coulomb EFG    0.47    (partially correlated — both sum
                                     point sources)
McConnell vs Ring Chi       0.41    (near-random — bonds vs rings)
Coulomb EFG vs Ring Chi     0.40    (near-random — different source
                                     geometries)
Random 5D vectors           0.36    (baseline)
BS vs HM                    0.999   (parallel — same physics, different
                                     mathematical approximation)
```

Different ring types' T2 patterns are NOT random pairs — they are
correlated because they share atomic positions. But the fitted
intensities differ (Phe -12 vs Tyr -11.28 vs Trp benzene -12.48 vs
Trp pyrrole -6.72), so the *parameterised* shielding contribution
per ring type is distinguishable in the calibration.

**Investigation implication.** The per-calculator report's "Pro
pyrrolidine impact analysis" must articulate, with the calculator's
specific physics, why Pro Ring with `Intensity = 0` produces zero
T2 contribution (linearity of the kernel × intensity product), AND
under what circumstances the calculator might emit a non-zero T2
value (any geometric kernel evaluated at atoms near Pro will produce
a non-zero raw kernel; only after intensity-weighting does the
shielding contribution zero out). Distinguish raw kernel from
shielding contribution.

### 4. Atom cyclic order determines ring-normal sign convention

From `Ring::ComputeGeometry` (`src/Ring.cpp`) and the sign-convention
verification in `PATTERNS.md` §22:

The ring normal is computed via SVD of vertex positions. The SVD's
right-singular vector for the smallest singular value is the
normal direction (up to sign). Sign is fixed by cross-product
convention: `n̂ = normalize((p₁ - p₀) × (p₂ - p₀))` where `p₀, p₁, p₂`
are the first three vertex positions in `atom_indices` order.

This means: reversing `atom_indices` reverses the cyclic walk
direction, which reverses the cross product, which flips `n̂`.

The shielding sign convention from `σ_ab = -dB_a^sec / dB_{0,b}` is
absorbed into:

```text
G_ab (Biot-Savart)    = -n_b · B_a · PPM_FACTOR
G_ab (Haigh-Mallion)  = -n_b · (H · n)_a
```

With these conventions, σ = I × G gives the correct physical sign:
literature I < 0 (diamagnetic) plus a probe atom above the ring
plane → σ > 0 (shielded, lower observed chemical shift).

The convention was VERIFIED analytically: Phe ring with I = -12.0,
probe 3 Å axial → σ = +1.40 ppm shielded (matches Case 1995). The
sign was wrong in early development; the bug was caught only by
analytical test, not by compilation or unit tests.

**Investigation implication.** Bundle C MUST preserve `atom_indices`
cyclic ordering bit-identically for non-Pro residues. Pro Ring's
ordering can be defined locally (no current calculator consumes Pro
Ring atoms for normal direction; future calculators define their
own conventions). For PHE/TYR/TRP-6 specifically, the substrate's
typed `RingPositionLabel` walk-order (`Ipso → Ortho1 → Meta1 → Para
→ Meta2 → Ortho2`) must match today's cyclic order produced by
`AminoAcidType::rings[].atom_names`. The investigation must answer:
does the substrate's `Ortho1` correspond to today's `CD1` or `CD2`?
Both possibilities exist depending on the substrate's labelling
convention; they are the same cyclic walk in opposite directions and
produce opposite-sign normals. Bundle C's verification gate is
shielding NPY bit-identity on bless fixtures; if the substrate's
labels are inversely-handed from today's atom_names, ConstructRings
FromSubstrate's `OrderedAtomIndices` helper must reverse to match.

### 5. HIS tautomer variants differ chemically

From `Ring.h` HID/HIE/HIP subclass docstrings + `topology-encoding-
dependencies-2026-05-05.md` §C.4 + literature:

Three protonation tautomers of histidine occur in proteins:

- **HID** (Nδ1-protonated): H on Nδ1; Nε2 is the bare-N. Neutral.
  Common in alkaline environments and as H-bond donor on the δ side.
- **HIE** (Nε2-protonated): H on Nε2; Nδ1 is the bare-N. Neutral.
  AMBER ff14SB default for "HIS" (per dependencies §D.1) because
  Nε2 protonation is more common in surveyed PDB structures.
- **HIP** (imidazolium, both NH protonated): formal charge +1 on
  Nε2 (per §C.4 — Lewis-localised; Nδ1 has formal charge 0; ring
  is symmetrised at the electron-density level but typed-localised
  at the formal-charge level). Common at acidic pH or when HIS is
  in a salt bridge.

Chemical-shift consequences (from Bertini, Luchinat, Parigi 2001 and
subsequent literature):

- HID Hε1 ~7.4 ppm; HIE Hδ2 ~7.2 ppm; HIP both Hδ2 + Hε1 ~8.5-9.0
  ppm (charged ring shifts ring protons strongly downfield).
- HID Cδ2 ~120 ppm; HIE Cδ2 ~135 ppm (15 ppm difference based on
  which N is protonated). This is a tautomer-distinguishing
  shielding signal; classical calculators must differentiate.
- H-bond donor/acceptor topology: HID donates from Nδ1, HIP donates
  from both, HIE donates from Nε2 only. NMR-assignment software
  uses this to distinguish.

In this project: ring-current intensity I = -5.16 nA/T is the
SAME across the three variants (per `Ring.h`); the variants differ
in which ring N is the H-bond donor, which N is bare, and the
imidazolium symmetry case for HIP. Ring chemistry shifts come not
from intensity differences but from per-atom RingPositionLabel
distinctions:

- HID: Nδ1 = `Heteroatom_NH`; Nε2 = `Heteroatom_NoH`.
- HIE: Nδ1 = `Heteroatom_NoH`; Nε2 = `Heteroatom_NH`.
- HIP: both `Heteroatom_NH`.

These per-position labels are the substrate's typed encoding of the
tautomer-distinguishing chemistry that downstream calculators must
respect.

**Investigation implication.** Each per-calculator report verifies
that the calculator does not pool HID/HIE/HIP atoms — the ring-
type-specific subclass dispatch via `Ring::TypeIndexAsInt()` should
distinguish them, but the agent must check the calculator code to
confirm. If a calculator pools HIS variants (e.g. uses
`HisImidazoleRing` indiscriminately), that's a finding for Bundle C.

### 6. TRP perimeter is real chemistry

From PATTERNS.md §24 and Case (1995):

Tryptophan's indole sidechain is a fused 5-ring (pyrrole) and
6-ring (benzene) system. The conjugated π electron current flows
around all 9 atoms of the indole perimeter, NOT independently in
each sub-ring. Modeling TRP as just `TrpPyrrole` + `TrpBenzene`
undercounts: the perimeter contribution is what the conjugated π
does as a single circuit.

Per-ring intensities sum exactly:

```text
I(TrpPerimeter) = I(TrpPyrrole) + I(TrpBenzene) = -6.72 + (-12.48)
                = -19.2 nA/T
```

This linearity is verified empirically (PATTERNS.md §24): batch
validation across 465 proteins shows BS-perimeter ratio = 1.000 and
HM-perimeter ratio = 1.000 after the `RingBondedExclusionFilter`
excludes shared-edge atoms from through-space evaluation. (The
13% HM excess at ratio 1.127 reported in earlier work was an
artifact from evaluating the surface integral at ring vertices
inside the source distribution; corrected 2026-04-02.)

The substrate currently encodes only `Indole_Trp_5` and
`Indole_Trp_6`; the 9-atom perimeter is a runtime construction in
`DetectAromaticRings` (`AminoAcidType.cpp:188` literal list:
`CG → CD1 → NE1 → CE2 → CZ2 → CH2 → CZ3 → CE3 → CD2`).

**Investigation implication.** The investigation must enumerate
every consumer of `RingTypeIndex::TrpPerimeter` and articulate how
Bundle C reproduces the perimeter Ring object. Two approaches
(documented in Bundle C planning):
- Synthesize at runtime in `ConstructRingsFromSubstrate` (no
  substrate change; perimeter atom-list union with defined cyclic
  order).
- Extend substrate with `RingSystemKind::Indole_Trp_9` (substrate-
  level encoding; generator must emit perimeter atoms with secondary
  RingPosition).
The investigation does not propose between approaches; it documents
what each consumer would need under each.

### 7. Pro pyrrolidine is silent today but real chemistry

From Vega & Boyer (1979) "Conformational analysis of proline residues
in a peptide chain," Biopolymers 18, 1797-1809; Sarkar et al. (1986)
"Conformational dependence of NMR chemical shifts in proline"; and
modern reviews (Schubert, Buynak, Schweitzer-Stenner 2002):

Proline's pyrrolidine ring puckers in two conformations:

- **C-γ exo**: Cγ above the Cα-Cδ-N plane.
- **C-γ endo**: Cγ below.

The endo/exo equilibrium shifts under different environments
(crystal vs solution; folded vs unfolded; coil vs β-turn). The
puckering is the residue's primary degree of freedom for chemical-
shift modulation:

- Cβ chemical-shift difference between exo and endo: ~1-2 ppm (¹³C).
- Cα chemical-shift difference: 0.3-0.8 ppm (¹³C).
- Cγ shows the largest exo/endo difference.
- Cα and Cβ shifts in trans-Pro vs cis-Pro residues additionally
  differ by 2-4 ppm (cis/trans isomerism around the preceding
  peptide bond) — a ring-puckering interaction with the upstream
  φ angle.

Currently no calculator in the project produces a per-frame Pro
puckering classification; no Pro chemistry signal flows through
calculator output beyond what backbone phi/psi already capture.
But the substrate has `RingSystemKind::Pyrrolidine_Pro` populated
for Pro N, Cα, Cβ, Cγ, Cδ. Materialising as a `Ring` object
exposes:

- `Ring::ComputeGeometry(positions)` returns ring centre, normal,
  vertices — natural inputs for puckering descriptors.
- Future trajectory-scope calculators can compute endo/exo
  classification via dihedral analysis of Cα-Cβ-Cγ-Cδ.
- Future ring-flip detection at Pro ring inversion (rare, ~kT
  barrier, occurs in unfolded states) becomes a per-frame query
  on the Pro Ring's geometry.

If Bundle C does NOT materialise Pro as a Ring object, the substrate
chemistry is structurally inaccessible at runtime — the
typed-substrate work loses concrete payoff for Pro. The user's
"no blur" rule explicitly addresses this: Pro Ring is materialised
NOT because current calculators consume it but because the substrate
encodes its existence as chemistry; the runtime must reflect that.

Reference: Joule, J. A. & Mills, K. (2010) "Heterocyclic Chemistry,"
5th ed., Wiley, ch. 7 — saturated heterocycles do not carry
circulating π current (intensity = 0 by physics), but their geometric
ring properties (ring strain, puckering modes) are real chemistry.

**Investigation implication.** Each per-calculator report's "Pro
pyrrolidine impact analysis" answers: with `Intensity = 0`, what
does the calculator produce per Pro Ring? Most ring-current
calculators: zero contribution to shielding via I × G. But the
calculator may produce GeometryChoice records (Pro ring iterated,
zero contribution recorded) — that's documentary value; not a
defect. Some calculators (DispersionResult uses ring atoms as
dispersion sources, not as ring-current sources; the dispersion
contribution is independent of intensity) may produce non-zero
contributions. The agent must trace each calculator's actual
behaviour from code, not assume intensity=0 means zero output.

### 8. Per-position chemical shifts differentiate at the substrate's granularity

From standard ¹³C NMR reference data (BMRB; Wuethrich 1986;
Cavanagh, Fairbrother, Palmer, Skelton 2007 "Protein NMR
Spectroscopy" 2nd ed.):

```text
Phe ring carbons (¹³C, ppm vs DSS):
  Cδ (ortho1, ortho2)       128-131
  Cε (meta1, meta2)          128-130
  Cζ (para)                  128-129
  Cγ (ipso)                  138-139

Tyr ring carbons:
  Cδ (ortho1, ortho2)        132-134
  Cε (meta1, meta2)          117-119  ← ortho-OH effect on meta C
  Cζ (para, with OH)         155-158  ← strong ipso-OH effect
  Cγ (ipso)                  127-130

Trp benzene carbons:
  Cε3 (ortho1)               118-120
  Cζ2 (ortho2)               113-115
  Cζ3 (meta1)                121-122
  Cη2 (meta2)                124-125

Trp pyrrole carbons:
  Cδ1 (PyrroleBeta)          124-126
  Cε2 (BridgeFusion)         137-139
  Cδ2 (BridgeFusion)         127-128  ← double bridgehead
```

The 5-15 ppm variation across ring positions WITHIN a single
residue is the angular fingerprint that the thesis's T2 calibration
must reproduce. Pooling positions ("aromatic Tyr carbon mean") loses
the para-OH effect at Cζ (157 ppm — strongly perturbed by OH) vs
Cε (118 ppm — less perturbed). The substrate's
`RingPositionLabel::Para` vs `Meta1` vs `Ortho1` typed labels are
the project's codification of this distinction.

**Investigation implication.** Each calculator's per-Ring outputs
that aggregate across ring positions (e.g. summed contributions over
all ring atoms without per-position decomposition) lose this
information. The investigation flags any such pooling as a finding
for Bundle C consideration; if a calculator currently pools, Bundle
C does not change the calculator (out of scope) but documents the
opportunity for future per-position-stratified calculator extension
once the substrate is the runtime authority.

### 9. Calibration depends on reproducible ring identity

From `learn/CLAUDE.md` and `spec/MATHS_GOALS.md`:

The thesis calibration pipeline (Python ridge regression on 720
proteins, 446K atoms, 55 kernels) reads the project's NPY output
arrays and fits per-element-and-atom-type stratified parameters.
R² = 0.818 settled 2026-04-10. The calibration's per-ring features
are keyed by `RingTypeIndex` enum values, which compile-time-pin to
the ring chemistry. Adding `RingTypeIndex::ProPyrrolidine` extends
the feature vector by one column (Pro ring's per-atom kernel
contribution). With Pro `Intensity = 0`, that column is identically
zero across all atoms — calibration sees zero variance in that
feature and fits a zero coefficient.

Adding the Pro column is not pollution: it preserves the
substrate-thesis correspondence that says "every ring chemistry
in the substrate gets a feature column." Future calculators that
extract puckering or ring-flip information add new columns
(`pro_pucker_endo_fraction`, `pro_ring_normal_dot_chi2`) that the
calibration can fit non-trivially. The Pro Ring object is the
addressing convention — if the Ring doesn't exist, the columns
can't be sized or named.

**Investigation implication.** The investigation report includes
the calculator's NPY output column count today, and what it would
become if Pro Ring is added. Most ring-current calculators emit
per-type columns sized by `RingTypeIndex::Count`. If `Count` is
hardcoded as `8`, Bundle C extends to `9` (or uses the enum
constant); the inventory must enumerate hardcoded sites.

### 10. GeometryChoice records make per-atom calculations reproducible

From `spec/GEOMETRY_CHOICE_BRIEF.md` (every agent reads this in
full):

A GeometryChoice is a runtime record of one decision: the calculator,
the entities involved, the role each played, the named numbers
(distances, intensities, thresholds) with units, and the outcome
(included / excluded / triggered with reason). Created during
`Compute()` via `GeometryChoiceBuilder::Record(...)` lambdas.

Concrete example from `BiotSavartResult`'s expected pattern (the
investigation will quote actual Builder calls verbatim):

```text
Choice #4321 — Biot-Savart at Lys33-HZ2 from Phe45 ring
  entities:
    field_atom    = Lys33 HZ2 (Element=H, Locant=Zeta, BranchAddress={1,0})
    source_ring   = Phe45 (RingTypeIndex=PheBenzene)
  named numbers:
    distance      = 6.4 Å
    intensity     = -12.0 nA/T
    near_field_threshold = 1.4 Å (from DipolarNearFieldFilter)
  outcome:
    included
  contribution:
    σ_T0 = +0.34 ppm; T2 = (-0.02, +0.18, -0.07, +0.04, -0.01) ppm
```

The methods chapter of the thesis cites these records as the
audit trail. Every shielding number in every calibration plot
traces back to its constituent decisions: which rings were
considered, which atoms were source vs probe, which filters
rejected which evaluations. Without records, the calibrated
shielding is a black box; with records, every prediction is
inspectable, reproducible, falsifiable.

When Pro rings join the system, calculators iterate them and
record outcomes (most likely "Pro ring iterated; Aromaticity = None;
no ring-current contribution") — that record is part of the audit
trail, not pollution. Future readers asking "why is the predicted
shift here +X ppm?" can answer "Pro residue Y was at distance Z;
Aromaticity = None excluded it from ring-current contribution per
filter Q." The audit gains coverage; the prediction is unchanged.

**Investigation implication.** Each per-calculator report quotes
the calculator's actual `Builder::Record(...)` lambda body verbatim,
identifies what gets recorded per ring iteration, and predicts what
records the calculator will emit for Pro rings. If the calculator
emits NO records (i.e. doesn't use `GeometryChoiceBuilder`), that's a
finding — the calculator violates PATTERNS.md "When implementing a
calculator" discipline.

### 11. Fail-loud is the scientific discipline, not an engineering preference

From this project's `feedback_no_attach_lifecycle_for_invariant_data`,
`feedback_no_summaries_as_gospel`, and the architecture's
substrate-miss FATAL pattern:

Thesis claims rest on completeness. A silent skip means a residue's
contribution is missing from a calculation; a future reader (a
reviewer, a replication group, a downstream calculator builder) may
not notice. A FATAL forces investigation: either fix the substrate
(extend chemistry coverage), document the load path's exclusion
(scope claim), or refuse the input (not all inputs are supported).

Bundle B already established this discipline: `ComposeAtomSemantic`
aborts on substrate misses; `RecanonicaliseAfterProtonation`'s
absence-of-rule means the validator catches non-canonical output.
Bundle C inherits the discipline: substrate-driven ring construction
on a load path that produces unrecognised ring chemistry must abort,
not silently degrade.

**Investigation implication.** Each per-calculator report identifies
the calculator's behaviour when ring chemistry is unrecognised
(e.g. a residue with `RingSystemKind::NotInRing` for atoms the
calculator expected to be in a ring). If the calculator fails-loud
(via assertion, FATAL+abort, or empty result + diagnostic), good.
If it silently produces a default, that's a finding for Bundle C —
substrate-driven construction should not silently introduce a
weaker invariant than today's string-matched code.

## Planned-calculator discipline — leave no chemistry data on the floor

Bundle C cannot be designed only for current ring-touching calculators.
The planned-calculator catalog (`spec/PLANNED_CALCULATORS_2026-04-22.md`,
`spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md`,
`spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md`) plus the just-
completed substrate audit at
`spec/plan/planned-calculator-substrate-audit-2026-05-06.md` enumerates
~60 planned calculators. Many consume ring chemistry — directly via
`Ring::*` virtuals + `Ring::atom_indices`, indirectly via per-atom
`ring_position` substrate fields, or via per-ring per-frame
trajectory-scope aggregates not yet built.

The user's discipline (memory entry `feedback_planned_calculators_stay_planned`):
**all proposed calculators are kept on principle**. Triage passes never
delete the queue. Bundle C must therefore not narrow the ring-data
surface such that planned calculators have to re-derive what they
need. Where the choice is between "narrow surface that satisfies
current calculators" and "richer surface that current calculators
ignore but planned calculators consume," default to the richer
surface unless there is a concrete cost.

### Concrete constraints this places on Bundle C design

The investigation reports surface these in their Q10 (Recommendations
for Bundle C); the framing here is the design lens.

#### A. Pro Ring atom ordering must support puckering descriptors

Pro puckering classification (endo vs exo at Cγ) needs the dihedral
N-Cα-Cβ-Cγ-Cδ. A future `ProPuckeringResult` ConformationResult
reads `Ring::atom_indices` and computes the dihedral.

**Constraint:** Bundle C's Pro Ring atom_indices walk is `N → Cα →
Cβ → Cγ → Cδ`, in residue-walk order. Not `Cα → N → ...` or any
other rotation. The cyclic walk's start atom + direction must be
stable across all proteins for downstream calculators to compute
the dihedral consistently. Documented in the ProPyrrolidineRing
class header.

Reference for puckering chemistry: Vega & Boyer (1979) Biopolymers
18, 1797-1809; Schubert, Buynak, Schweitzer-Stenner (2002) for
modern classification methodology.

#### B. RingPositionLabel must be reachable per atom from runtime code

Per the substrate audit (Section 2 / Section 5), planned per-position
stratification calculators consume `RingPositionLabel` per atom:

- **Per-SS CSA stratification diagnostic** (PLANNED_CALCULATORS §3):
  stratifies CSA tensors by (backbone_role, ring_position) cross.
  Reads `LegacyAmber().SemanticAt(ai).ring_position.primary.position`
  per ring atom.
- **PelloniDifferentialBiotSavart** and **DistributedRingCurrent**
  (NMR_EXTRACT_DESIDERATA §B) compute per-position contributions to
  shielding instead of per-ring totals; need `RingPositionLabel`
  per ring vertex.
- **Per-position aromatic ¹³C reference** diagnostics: Tyr Cζ vs Cε
  vs Cδ stratification (40 ppm spread per ring per the science
  section above) needs reachable per-position labels.

**Constraint:** Bundle C does not block this. Substrate already
exposes `ring_position` per atom via `SemanticAt(ai)`. The Ring
object itself does NOT need a per-vertex `RingPositionLabel` member
(redundant with substrate). Calculators that want per-position info
read substrate via atom_indices. Bundle C confirms this access path
works for all 9 ring types after substrate-driven construction;
agents check that `SemanticAt(atom_indices[k]).ring_position.primary.position`
returns the expected label for each k for fixture residues.

#### C. TRP perimeter representation choice has planned-calculator consequences

Two options for Bundle C (named earlier):
- (b) Synthesize at runtime in `ConstructRingsFromSubstrate`. No
  substrate change.
- (a) Extend substrate with `RingSystemKind::Indole_Trp_9`.

Planned-calculator implications:

- **Per-ring trajectory aggregates** (RingPlanarityWelford,
  RingNormalDriftWelford — NMR_EXTRACT_DESIDERATA §A.10): work
  identically under either approach because they iterate
  `protein.RingAt(ri)` and the perimeter Ring exists in either case.
- **Per-atom-per-ring trajectory stats** (RingNeighbourhoodTrajectoryStats):
  same — atom_indices is the input regardless.
- **Per-position-on-perimeter calculators** (hypothetical: "shielding
  contribution from atoms at Trp perimeter ortho1 vs meta2 across
  the trajectory"): under (a), the substrate emits perimeter atoms
  with secondary `RingPosition.position` tagged with perimeter labels;
  under (b), runtime calculators have to derive these from the
  union of TrpPyrrole + TrpBenzene atoms with the cyclic walk
  defined locally. (a) is more reusable.
- **Substrate-driven ring construction itself**: under (a), perimeter
  emerges naturally from the per-atom `ring_position.secondary` field
  (atoms in two rings emit two records); under (b), `ConstructRingsFromSubstrate`
  has perimeter synthesis logic embedded.

Bundle C should evaluate both with this lens. The investigation does
not propose between them; it surfaces the planned-calculator angle
so Bundle C drafting weighs both consequences.

#### D. Ring object stability across frames matters for trajectory-scope calculators

Trajectory-scope calculators (`PerFrameExtractionSet` and successors)
iterate rings per frame. The `Ring::atom_indices` and
`Ring::parent_residue_index` are invariant; only `Ring::accumulated`
mutates. Bundle C must not introduce per-frame ring re-construction
or ring-identity churn — the rings_ vector lives on Protein (the
invariant) and is built once at `Protein::FinalizeConstruction`.

**Constraint:** Bundle C's `ConstructRingsFromSubstrate` runs once
at `FinalizeConstruction`. Pro Ring objects are built at the same
time as aromatic Rings. Trajectory-scope calculators see the same
ring inventory at every frame.

This is already the architectural guarantee for current calculators;
flagged here so the investigation confirms no current calculator
re-builds rings per frame (would be a bug today).

#### E. Calibration feature-vector schema is sized by RingTypeIndex::Count

Per the science section §9: per-element-and-atom-type ridge regression
on 55 kernels at R² = 0.818. Per-ring-type features key on
`RingTypeIndex` enum integer values. Adding `RingTypeIndex::ProPyrrolidine`
extends the feature vector by one column.

**Constraint:** Bundle C uses `static_cast<size_t>(RingTypeIndex::Count)`
or `kRingTypeCount` for array sizing wherever ring-type-keyed arrays
appear (`per_type_G_T0_sum[Count]` on ConformationAtom et al.). If
`Count` is hardcoded as `8` anywhere, Bundle C extends to `9`. The
investigation enumerates all hardcoded `[8]` sites in a calculator-
specific report (Q6 of the per-report template).

Planned-calculator implication: calibration code (`learn/c_equivariant/`)
needs to know that `RingTypeIndex::ProPyrrolidine`'s column is
identically zero. The downstream calibration pipeline reads
project's NPY files via `python/nmr_extract/_catalog.py`; that
catalog needs an entry for the Pro column once Bundle C lands.
Coordination with `python/` SDK team is part of Bundle C's
deliverable, not Phase 2.

#### F. Hypothetical per-residue ring catalog for non-standard residues

The substrate generator currently covers the standard 20 + 10
variants. Planned work (mentioned in NMR_EXTRACT_DESIDERATA §F):
non-standard residues (modified amino acids — phospho-Ser, methyl-
Lys; nucleic-acid bases; cofactors). These have rings the current
substrate doesn't encode.

**Constraint:** Bundle C's substrate-driven ring construction must
gracefully handle residues whose substrate is empty (no `ring_position`
for any atom). The fail-loud invariant from §11 of the science
section says: if a residue is `AminoAcid::Unknown` with named
atoms, abort. But for non-standard residues that legitimately have
no current substrate coverage, the calculator should emit an empty
rings_ contribution for that residue and emit a GeometryChoice
record naming the unsupported case.

Planned-calculator implication: when non-standard substrate
extension lands (Phase 2 of substrate work), `ConstructRingsFromSubstrate`
extends naturally — new RingSystemKind values, new RingTypeIndex
values, new substrate variant tables. Bundle C should not bake
"ring chemistry is exclusively the standard 20" into its
construction logic; the construction is data-driven from substrate.

### Adjusted per-report template — Q10 expansion

Question 10 ("Recommendations for Bundle C") now has two parts:

**Q10.a — For current calculators.** Specific, typed answers in
the form: "Calculator does X today; under Bundle C's substrate-
driven construction with Pro Ring at Intensity=0, calculator does
Y; Bundle C must / should / need not Z." Concrete file:line
references.

**Q10.b — Implications for planned calculators that share this
calculator's ring chemistry surface.** Reference the planned-
calculator audit at `spec/plan/planned-calculator-substrate-audit-
2026-05-06.md`; identify which planned calculators consume this
calculator's outputs OR consume the same ring substrate fields OR
implement related physics; state how Bundle C decisions affect
their feasibility. Surface "data left on the floor" cases — where
narrowing the ring-data surface would block planned work.

### Reading list addition

Every agent additionally reads:

- `spec/plan/planned-calculator-substrate-audit-2026-05-06.md` —
  the just-completed mapping of planned calculators to substrate
  fields. Required reading for Q10.b.
- `spec/PLANNED_CALCULATORS_2026-04-22.md` — five Session-0 ideas.
- `spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md` — ~40 calculator/I/O/
  diagnostic ideas; sections A (new ConformationResults), B
  (variations on existing kernels), C (I/O surfaces), D (diagnostics).
- `spec/PLANNED_CALCULATORS_TIME_SERIES_2026-04-24.md` — paper→
  candidate-calculator entries.

Memory entry `feedback_planned_calculators_stay_planned` is auto-
loaded; it codifies the discipline.

## Typed-identity discipline (CRITICAL — read before code)

**Every atom reference in every investigation report is in typed
identity terms.** Strings appear ONLY as verbatim quotes from
existing code, immediately followed by their typed-identity
translation. Findings, recommendations, calculator-impact analysis
— all in typed vocabulary. The substrate's typed enum vocabulary IS
the noun set for reports.

The substrate vocabulary, in `src/SemanticEnums.h`:

```text
Element              :: H, C, N, O, S, Unknown
Locant               :: None, Alpha, Beta, Gamma, Delta, Epsilon, Zeta, Eta
BranchAddress        :: { outer ∈ 0..2, inner ∈ 0..2 } 2-level disambiguator
DiastereotopicIndex  :: None, Position2, Position3
BackboneRole         :: None, Nitrogen, AlphaCarbon, CarbonylCarbon,
                        CarbonylOxygen, AmideHydrogen, AlphaHydrogen
TerminalState        :: Internal, NtermCharged, NtermNeutral,
                        CtermDeprotonated, CtermProtonated
RingSystemKind       :: NotInRing, Benzene_Phe, Benzene_Tyr,
                        Imidazole_His, Indole_Trp_5, Indole_Trp_6,
                        Pyrrolidine_Pro
RingPositionLabel    :: NotInRing, Ipso, Ortho1, Ortho2, Meta1, Meta2,
                        Para, PyrroleAlpha, PyrroleBeta, BridgeFusion,
                        Heteroatom_NH, Heteroatom_NoH, Heteroatom_OH,
                        Saturated
RingMembership       :: { ring: RingSystemKind, position:
                          RingPositionLabel, ring_size, aromatic,
                          planar, n_heteroatoms }
RingPosition         :: { primary: RingMembership, secondary:
                          RingMembership } (secondary populated for
                          fused atoms)
PolarHKind           :: NotPolar, BackboneAmide, SidechainPrimaryAmide,
                        IndoleNH, AmmoniumNH, GuanidiniumNH,
                        ImidazoleNH, CarboxylOH, HydroxylOH_Aliphatic,
                        HydroxylOH_Aromatic, ThiolSH, AmineNH,
                        OtherPolarH
PlanarGroupKind      :: None, PeptideAmide, SidechainAmide,
                        Guanidinium, Imidazole, Aromatic6Ring,
                        Aromatic5Ring, Carboxylate, AromaticHydroxyl,
                        AromaticOxide
ProchiralStereo      :: NotProchiral, ProR, ProS, Unassigned
PlanarStereo         :: NotApplicable, E, Z, Unspecified
PseudoatomKind       :: None, M, Q, R
AtomMechanicalIdentity :: { Element, Locant, BranchAddress,
                            DiastereotopicIndex, BackboneRole }
AtomSemanticTable    :: identity 5-tuple plus 9 chemistry fields
                        (prochiral, planar_group, planar_stereo,
                        pseudoatom, polar_h, ring_position, aromatic,
                        formal_charge, is_exchangeable,
                        equivalence_class)
```

When the existing code reads `atom.pdb_atom_name == "CG"` and
intersects with `AminoAcidType::rings[].atom_names`, the report
describes the atom as **"the ipso carbon of the Phe ring (Element=C,
Locant=Gamma, BranchAddress={0,0}, RingPositionLabel=Ipso,
RingSystemKind=Benzene_Phe)"** — citing the verbatim string from the
code, then immediately translating. The agent's recommendations and
calculator-impact analysis use the typed terms only; the strings
exist as historical artifacts being deleted.

If an investigation report puts strings in its FINDINGS or
RECOMMENDATIONS sections (anywhere outside verbatim code quotes),
that is a documentary regression and must be flagged + rewritten.

**No shortcut framing. No "this atom is HD1" without "(Element=H,
Locant=Delta, BranchAddress={1,0}, in Imidazole_His)".**

## String-discipline check (a mandatory finding in every report)

Every report includes a section titled "String discipline check"
that names every direct read of `atom.pdb_atom_name` or substring
operations on atom names that the calculator performs. If the
calculator uses `AtomMechanicalIdentity` typed lookups, say so. If
any string read appears, that's a finding for Bundle C to clear —
flag it loud.

The investigation does NOT require the agent to fix any string
reads; it requires the agent to enumerate them so Bundle C planning
addresses them.

## Science context — ring chemistry

The investigation report needs to describe each calculator's
physics, not just its code shape. Eight ring types currently exist
in the project:

```text
PheBenzene     :: 6-ring aromatic. I = -12.0 nA/T (Giessner-Prettre 1969).
TyrPhenol      :: 6-ring aromatic with para-OH. I = -11.28 nA/T.
TrpBenzene     :: 6-ring of the indole (fused with pyrrole). I = -12.48 nA/T.
TrpPyrrole     :: 5-ring of the indole, fused with benzene. I = -6.72 nA/T.
TrpPerimeter   :: 9-atom synthetic perimeter ring (the conjugated π
                  current encircling the indole). I = -19.2 nA/T (sum of
                  TrpPyrrole + TrpBenzene). Built at runtime as the
                  union of TrpPyrrole + TrpBenzene atoms minus the
                  shared edge, with a defined cyclic ordering.
HisImidazole   :: 5-ring with two nitrogens. Used when protonation state
                  is ambiguous (variant_index = -1). I = -5.16 nA/T.
HidImidazoleRing :: HIS variant 0 (Nδ1-protonated tautomer).
HieImidazoleRing :: HIS variant 1 (Nε2-protonated tautomer; AMBER ff14SB
                    default for HIS, per dependencies §D.1).
```

Bundle C adds:

```text
ProPyrrolidine :: 5-ring saturated heterocycle (N-Cα-Cβ-Cγ-Cδ). I = 0
                  (no ring current, saturated). Cited per Joule & Mills
                  "Heterocyclic Chemistry" 5e (2010) ch. 7 — saturated
                  heterocycles do not carry circulating π current.
                  Aromaticity = None, NitrogenCount = 1, RingSize = 5.
```

Per-calculator physics (each agent describes their own; this list is
orientation):

- **Biot-Savart** (`BiotSavartResult`): Johnson-Bovey two-loop line
  integral. B-field at probe atom = ∫ (I dℓ × r̂) / r² over the ring
  vertices. Each ring vertex contributes a line segment. Sign
  convention: `G_ab = -n_b · B_a · PPM_FACTOR` so that with
  literature intensities (I < 0 diamagnetic), σ = I·G gives the
  correct shielded sign above the ring plane. Ring-normal direction
  matters; SVD on vertex positions; ordering of `atom_indices`
  determines normal sign via cross-product convention.

- **Haigh-Mallion** (`HaighMallionResult`): surface-integral over
  ring face, evaluated via 7-point Gaussian quadrature (Stroud
  T2:5-1) with adaptive subdivision near the ring face (level 1 at
  2.0 Å; level 2 at 1.0 Å; max 2 levels). Same sign convention via
  `G = -n⊗V` where V = H·n.

- **Pople ring susceptibility** (`RingSusceptibilityResult`):
  point-dipole at ring center with anisotropic susceptibility tensor.
  σ ∝ (3 n̂_a n̂_b - δ_ab) χ_ring / r³. Per ring type's Δχ.

- **Pi-quadrupole** (`PiQuadrupoleResult`): ring's quadrupole moment
  at center + per-vertex distribution; 1/r⁹ leading term.

- **Dispersion** (`DispersionResult`): London dispersion through ring
  atoms with CHARMM switching function (Brooks 1983) tapered between
  R_switch = 4.3 Å and R_cut = 5.0 Å.

The investigation report cites each calculator's physics with
references; the substrate-Pro-Ring impact is then derivable: with
`Intensity = 0`, every ring-current calculator produces zero
contribution from Pro; with `Aromaticity = None`, filters that gate on
aromatic should reject; with `RingSize = 5` matching pyrrolidine,
distance/extent computations remain meaningful even though
contribution is zero.

## Per-report template

Each per-calculator report at
`spec/plan/ring-investigation-2026-05-06/<CalculatorName>.md`
answers, in order:

### 1. What the calculator does (physics-first, not code-first)

The science of the calculation — what physical quantity is computed,
what equation, what units, what literature it cites. One paragraph
of prose, then equations in display form.

### 2. Ring-information consumption (typed)

Which `Ring::` virtual methods does the calculator's `Compute()`
read? Per-ring-type values it consumes (e.g. `ring.Intensity()` for
each `RingTypeIndex`). Whether it reads `ring.atom_indices` and
how it uses the indices (per-atom iteration; centroid; SVD-normal
input). Whether it reads `RingNeighbourhood` records. All atom
references in typed identity.

### 3. Iteration and filter shape

Per-atom-per-ring? Per-ring-only? Per-atom-with-ring-lookup? What
KernelFilterSet members does it apply? Each filter's physics
rejection criterion (Why does this filter exist? What kernel
breakdown does it prevent?). The exact rejection condition in
typed terms (e.g. "filter rejects when probe atom is one of the ring
vertices" rather than "filter rejects atoms in the same ring").

### 4. GeometryChoice records emitted

What does the calculator's `GeometryChoiceBuilder::Record(...)`
emit per ring iteration? Quote the lambda body verbatim. What
entities are attached (atom, ring, bond), with what roles
("ring_source", "field_point", etc.), what named numbers (distance,
intensity, threshold), what units? This is the audit-trail surface
Bundle C must understand.

### 5. Output surface

NPY columns emitted via `WriteFeatures` — list each one with shape +
dtype + chemistry meaning. ConformationAtom fields written — list
each with type and unit. Per-Ring accumulated state — list each
field on `Ring::accumulated`. H5 columns if applicable.

### 6. RingTypeIndex interaction

Every site using `RingTypeIndex::Count`, `[8]` array sized by ring
type, or `RingTypeIndex::*` enum values directly. Exact line
references. If hardcoded `[8]`, flag it for Bundle C's extension
work.

### 7. AminoAcidType::rings interaction

Whether the calculator reads `aatype.rings[]` directly or
`aatype.rings[].atom_names`. Bundle C deletes the latter; the
investigation confirms only `DetectAromaticRings` reads it (likely
true). Flag any other reader as a Bundle C migration target.

### 8. Pro pyrrolidine impact analysis

Given a Pro Ring with `Aromaticity = None`, `Intensity = 0`,
`RingSize = 5`, `NitrogenCount = 1`, what does this calculator's
`Compute()` produce when iterating it?
- Filter outcomes: which filters reject, which let through.
- Per-iteration contribution: zero (intensity-weighted) or whatever
  the math gives.
- GeometryChoice records emitted: typed quotes of the records that
  would land in the choice log.
- Calculator output surface impact: do any NPY fields,
  ConformationAtom fields, or per-Ring accumulators receive non-zero
  values? Should they?
The agent translates the typed inputs (substrate-derived) through
the calculator's physics to predict output behaviour — a real
physics description, not "may have effect."

### 9. String discipline check

Every direct read of `atom.pdb_atom_name`, `residue.three_letter`
substring, or other string operation in the calculator. List
verbatim with line references. Translate each to the typed-identity
equivalent. If any reads exist, flag for Bundle C clearance.

### 10. Recommendations for Bundle C

Specific, typed answers to the question "does Bundle C need to do
anything for this calculator?". Possibilities:
- "Calculator transparently iterates new Pro rings and produces
  zero contribution. No Bundle C work required. GeometryChoice
  records will gain Pro entries naturally — documentary value."
- "Calculator has hardcoded `[8]` array at line N; Bundle C must
  extend to `[9]` or `RingTypeIndex::Count`."
- "Calculator's filter passes Pro rings through computation; Pro
  contribution is mathematically zero by intensity but adds CPU.
  Consider an explicit `Aromaticity::None` skip-filter for
  performance — though see PATTERNS.md anti-premature-optimization."
- "Calculator reads atom-name string at line N; clearance required
  before / during Bundle C."

## Investigation process

1. **Inventory phase** (this dispatch). One Explore agent does a
   case-insensitive grep + targeted read across all calculator code
   to identify every site that touches ring data. Returns a list
   with one-line summaries per site, distinguishing substantive
   consumers from incidental references. Lands at
   `spec/plan/ring-investigation-2026-05-06/inventory.md`.

2. **Verification.** The user (or this assistant) runs an independent
   grep to confirm the inventory's coverage didn't tucker out.
   Patches gaps if any.

3. **Per-calculator phase.** One Explore agent per substantive
   consumer in the verified inventory. Each reads its calculator,
   answers the 10-question template, lands at
   `spec/plan/ring-investigation-2026-05-06/<CalculatorName>.md`.

4. **Codex parallel review.** User dispatches codex per calculator
   report. Each codex run examines the report against the calculator
   code; returns findings.

5. **Bundle C drafting.** With the verified inventory + per-
   calculator reports + codex findings, Bundle C's brief is drafted
   with concrete answers (no speculation).

6. **Iterate** if a per-calculator round surfaces new questions
   needing follow-on investigation.

## Output locations

```text
spec/plan/ring-investigation-2026-05-06/
    README.md                                  ← this file (framing memo)
    inventory.md                               ← inventory phase output
    <CalculatorOrSurfaceName>.md               ← per-calculator reports
```

Predetermined paths so codex's per-calculator parallel review can
target each report directly.

## Reading list (every agent reads, in order)

1. **This README.md.**
2. `src/SemanticEnums.h` — substrate vocabulary; the typed
   identity tuple definitions; every ring-related enum.
3. `src/Ring.h` and `src/Ring.cpp` — the Ring class hierarchy:
   base class virtuals, subclass overrides, RingGeometry,
   RingAccumulated.
4. `src/generated/LegacyAmberSemanticTables.h` — `LookupBy` /
   `LookupCap` / `ApplyCapDelta` / `ParseAtomName` /
   `ComputeAtomMechanicalIdentity`.
5. `src/LegacyAmberTopology.h` — `SemanticAt`, `IdentityAt`,
   `ResidueAtomsWithIdentity`, `AtomWithRole`.
6. `src/Protein.cpp` — `DetectAromaticRings` (the function being
   replaced), `FinalizeConstruction` ordering.
7. `src/AminoAcidType.cpp` — the `rings[].atom_names` const-char-
   pointer arrays for understanding how today's code identifies
   ring atoms.
8. **`spec/GEOMETRIC_KERNEL_CATALOGUE.md`** — full physics
   derivations of each calculator. Required physics grounding.
9. **`spec/GEOMETRY_CHOICE_BRIEF.md`** — what a GeometryChoice IS,
   the Builder API, the per-calculator manifest. Required for the
   choice-records section of each report.
10. `spec/CONSTITUTION.md` — sign conventions and inviolable rules.
11. `OBJECT_MODEL.md` Atom + Ring + ConformationAtom + RingNeighbourhood
    sections.
12. `PATTERNS.md` esp. §"When implementing a calculator", §"T2
    Completeness", anti-patterns, the rule-application architecture
    WIP section.
13. `spec/plan/topology-encoding-dependencies-2026-05-05.md` §H
    (substrate composition rules including the cap-overlay field-
    level discipline and the deletion-variant deny overlay).
14. `spec/plan/naming-applicator-architecture-sketch-2026-05-06.md`
    (naming-application architecture as recently landed).
15. The calculator's own source files for the per-calculator phase.

## Anti-asks (for every agent)

- Do NOT propose Bundle C's design. Bundle C is drafted after
  investigation lands. The investigation produces empirical answers,
  not architectural proposals.
- Do NOT modify any code. This is read-only investigation.
- Do NOT use atom-name strings in findings or recommendations.
  Strings appear only in verbatim code quotes with immediate typed-
  identity translation.
- Do NOT shorten the report to fit a length limit. The user's
  directive is "hardcore full knowledge of the physics and use" — no
  size cap.
- Do NOT assume calculator behaviour without reading the source.
  Every claim about "this calculator does X" cites the line in the
  source; "may," "might," "probably" language is forbidden.
- Do NOT skip the science section. Each report grounds the
  calculator's behaviour in physics, not just code shape.

## What the inventory phase produces

The inventory dispatch (next; first investigation agent) produces
`inventory.md` with:

- Case-insensitive grep results across `src/` for "ring",
  "Ring", "RING" — every line with file + line number + brief
  context.
- Targeted reads of each calculator file (`*Result.cpp`) for ring
  interactions.
- Per-site classification:
  - SUBSTANTIVE consumer (uses ring data for physics — needs full
    per-calculator report).
  - INCIDENTAL reference (mentions rings in comments, includes
    headers, accesses peripherally — does not need full report).
  - INFRASTRUCTURE (Ring class itself, DetectAromaticRings,
    CovalentTopology::Resolve, KernelFilterSet ring-related members
    — needs an infrastructure report, possibly bundled).
- One-line summary per site: which Ring API it uses, what
  computation it performs, whether it has any string reads on atom
  names.
- A coverage diagnostic: "I read N files and grep returned M
  matches; if you re-run the grep yourself you should see the same
  M matches." So the user can verify the agent didn't tucker out.

The inventory is the input to per-calculator dispatch decisions:
which calculators get a full report, which can be folded into the
infrastructure report, which need no investigation.

## Notes for context recovery

If this conversation runs out of context and a future session opens
cold, this README is the recovery anchor for the ring investigation.
Memory entry pointer:
`memory/project_bundle_c_ring_investigation_20260506.md`. The
investigation's products live at `spec/plan/ring-investigation-2026-
05-06/`. Bundle C drafting is gated on the investigation completing
and being codex-reviewed.
