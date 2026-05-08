# Substrate-completeness audit for Bundle C ring construction (2026-05-07)

Read-only audit of the chemistry substrate (`AtomSemanticTable`,
`RingTypeIndex`, `RingPositionLabel`, generator output, residue
reference) against the requirements of substrate-driven ring
construction (Bundle C) and the planned-calculator pipeline. Every
atom reference uses typed identity; every chemistry-driven gap
recommendation cites primary chemistry literature or BMRB / IUPAC
convention. Substrate-side scope only — no calculator-side
prescriptions.

Locked decisions per
`memory/project_bundle_c_ring_investigation_20260506.md` and the
2026-05-07 evening session are taken as contract; the audit
identifies completeness gaps in the substrate as it stands at
HEAD.

Sources read for this audit:

- `src/SemanticEnums.h` (substrate enum declarations)
- `src/Ring.h`, `src/Ring.cpp` (Ring-class hierarchy + ComputeGeometry)
- `src/Types.h:175-199` (RingTypeIndex enum + RingTypeName mapping)
- `src/generated/LegacyAmberSemanticTables.cpp` (substrate generator
  output: per-residue, per-variant tables)
- `src/generated/LegacyAmberSemanticTables.h` (LookupBy / LookupCap /
  ApplyCapDelta / ParseAtomName)
- `src/LegacyAmberTopology.h`, `.cpp` (substrate-access surface
  including SemanticAt and IdentityAt)
- `src/AminoAcidType.h`, `.cpp` (rings[].atom_names — string surface
  to be removed)
- `src/Protein.cpp:548-646` (`DetectAromaticRings` migration target)
- `tools/topology/build_semantic_tables.cpp` (substrate generator
  source; existence verified, code not transcribed)
- `spec/plan/topology-residue-reference-2026-05-05.md` (canonical
  per-residue, per-atom reference)
- `spec/plan/chemistry-question-3-ring-topology-2026-05-05.md`
  (Path A architectural ground truth)
- `spec/plan/planned-calculator-substrate-audit-2026-05-06.md`
  (cross-reference — what planned calculators will consume from
  the ring substrate)
- `spec/plan/ring-investigation-2026-05-06/inventory.md` (current
  consumers; per-calculator `[8]` hardcoding noted)


## Executive summary

The substrate is **construction-ready** for the locked Bundle C
shape on PHE, TYR, HIS (HID/HIE/HIP), and TRP-5/TRP-6 — the per-atom
`RingPosition.primary.{ring, position}` fields are populated for
every aromatic ring atom in every variant in the generated tables,
and the typed-position labels match Markley 1998 §2.1.1 and the
chi₂-priority convention for ortho/meta disambiguation. PRO is
also populated at `Pyrrolidine_Pro/{Heteroatom_NoH | Saturated}`
across N, Cα, Cβ, Cγ, Cδ; the locked Pro 5-label extension
(`ProRingNitrogen` / `ProRingAlphaCarbon` / `ProRingBeta` /
`ProRingPuckerPivot` / `ProRingDelta`) is **additive — not
implemented yet**.

Five substrate-side gaps prevent Bundle C and the projection /
planned-calculator slices from reading purely typed answers. In
priority order, with Pro labels being a locked decision the audit
treats as a contract not a debate:

1. **`RingTypeIndex` lacks `ProPyrrolidine` and `Count` is hardcoded
   `8`.** Adding the Pro index extends the enum to 9 values; the
   per-aromatic-type calculator arrays continue to gate on
   `< 8` (locked, per Bundle C scope) — the substrate gains the
   index, calculators are unchanged.
2. **`Ring` class hierarchy lacks `ProPyrrolidineRing`.** The
   factory `CreateRing(RingTypeIndex)` and the `FiveMemberedRing`
   sister classes are the obvious place for a Pro-specific
   subclass; locked decision.
3. **Pro per-atom `RingPositionLabel` enrichment is missing five
   typed values** (`ProRingNitrogen`, `ProRingAlphaCarbon`,
   `ProRingBeta`, `ProRingPuckerPivot`, `ProRingDelta`). Existing
   `Saturated` (value 13) and `Heteroatom_NoH` (11) on Pro atoms
   today are coarse and lose the puckering-relevant chemistry.
4. **TRP perimeter is not encoded in the substrate.** The 9-atom
   perimeter is a runtime synthesis from the union of `Indole_Trp_5`
   atoms + `Indole_Trp_6` atoms minus shared bridgeheads. No
   `RingSystemKind::Indole_Trp_9` exists; per-atom secondary
   `RingPosition.secondary` carries only the bridgehead-fusion
   relation, not perimeter ortho/meta/para labels.
5. **Ring atom_indices ordering convention is undocumented.** It
   exists de facto in `AminoAcidType.cpp` (PHE: `CG → CD1 → CE1 →
   CZ → CE2 → CD2`, etc.) and is load-bearing for the SVD ring-normal
   sign convention (PATTERNS.md §22). Bundle C's
   `OrderedAtomIndices` helper reads from substrate
   `RingPositionLabel`; the substrate's de facto ordering
   convention (Ipso → Ortho1 → Meta1 → Para → Meta2 → Ortho2 for
   benzene rings) needs to be **documented in `Ring.h` per
   ring chemistry**, with the Pro convention `N → Cα → Cβ → Cγ
   → Cδ` documented as locked.

Three smaller gaps rated below the top five:

6. Coverage gaps in `RingMembership.aromatic` / `.planar` semantics
   for the HIP variant (formal-charge 0 on Nδ1, formal-charge +1 on
   Nε2 per the Lewis-localised convention; both N's tagged
   `Heteroatom_NH`, but this loses the imidazolium charge
   asymmetry information).
7. No typed predicate at the residue level for "this residue's atoms
   collectively form which `RingSystemKind`s?" (today derived by
   iterating per-atom). Adding one to `LegacyAmberTopology` is
   convenience.
8. No typed projection from `RingPositionLabel` to the canonical
   IUPAC / BMRB atom name. AMBER names ARE the de facto storage on
   `Atom.pdb_atom_name`; IUPAC ≡ BMRB ≡ AMBER for ring atoms in
   the standard 20 (modulo Markley 1998 protonation-tautomer
   nomenclature for HID/HIE/HIP), so projection is shallow but
   not yet expressed as a typed function.

The full picture: the substrate is the right shape; the coverage is
nearly complete on the aromatic side; the Pro extension is
straightforward additive work; the perimeter and ordering
conventions are the conceptually deepest gaps because they involve
synthesised chemistry (perimeter) and load-bearing physics
(SVD-normal sign convention via cyclic order).


## Current substrate state per ring chemistry

This section walks each of the seven `RingSystemKind` values
populated by the generator today (eight `RingTypeIndex` values, but
HID/HIE/HIP all map to one `RingSystemKind::Imidazole_His` — so the
substrate has fewer ring-system kinds than ring-type indices),
documenting what's encoded and what's missing.

### Phe — `RingSystemKind::Benzene_Phe`

**Per-atom substrate (from `kPheAtoms`,
`src/generated/LegacyAmberSemanticTables.cpp:298-319`):** every ring
atom in the 6-ring carries the typed position label.

| AMBER name | typed identity | RingPosition.primary | aromatic |
|---|---|---|---|
| `CG` | `(C, Gamma, {0,0}, None, None)` | `Benzene_Phe / Ipso / 6 / true / planar=true / 0 het` | true |
| `CD1` | `(C, Delta, {1,0}, None, None)` | `Benzene_Phe / Ortho1 / 6 / true / true / 0` | true |
| `CD2` | `(C, Delta, {2,0}, None, None)` | `Benzene_Phe / Ortho2 / 6 / true / true / 0` | true |
| `CE1` | `(C, Epsilon, {1,0}, None, None)` | `Benzene_Phe / Meta1 / 6 / true / true / 0` | true |
| `CE2` | `(C, Epsilon, {2,0}, None, None)` | `Benzene_Phe / Meta2 / 6 / true / true / 0` | true |
| `CZ` | `(C, Zeta, {0,0}, None, None)` | `Benzene_Phe / Para / 6 / true / true / 0` | true |
| `HD1` | `(H, Delta, {1,0}, None, None)` | `Benzene_Phe / Ortho1 / 6 / true / true / 0` | false |
| `HD2` | `(H, Delta, {2,0}, None, None)` | `Benzene_Phe / Ortho2 / 6 / true / true / 0` | false |
| `HE1` | `(H, Epsilon, {1,0}, None, None)` | `Benzene_Phe / Meta1 / 6 / true / true / 0` | false |
| `HE2` | `(H, Epsilon, {2,0}, None, None)` | `Benzene_Phe / Meta2 / 6 / true / true / 0` | false |
| `HZ` | `(H, Zeta, {0,0}, None, None)` | `Benzene_Phe / Para / 6 / true / true / 0` | false |

Ring-attached H atoms also carry the typed ring-position label of
their parent C (mirrored for downstream H-position-stratified
calculators per planned-audit §B4 lines 264-272).

**Cyclic walk in substrate `RingPositionLabel` enum order**:
`Ipso(1) → Ortho1(2) → Ortho2(3) → Meta1(4) → Meta2(5) → Para(6)`
— but this is **not** the natural cyclic walk around the ring.
The chemistry-correct walk is `Ipso → Ortho1 → Meta1 → Para →
Meta2 → Ortho2`, matching `AminoAcidType.cpp:147` which lists
`{"CG","CD1","CE1","CZ","CE2","CD2"}`. Bundle C's
`OrderedAtomIndices` must encode this canonical walk; the enum
integer order is a label-listing order, not a cyclic-traversal
order.

**Planar-group coverage (`F7 planar_group`).** Every ring atom
carries `PlanarGroupKind::Aromatic6Ring` (per
`SemanticEnums.h:268`). All six C and the five attached H atoms
are aromatic.

**Markley 1998 / IUPAC alignment.** Locant labels match Markley
2.1.1 (Greek-letter side-chain locants); `BranchAddress.outer`
distinguishes CD1 vs CD2 / CE1 vs CE2 per the chi2-priority rule
(Markley 1998 Fig 1 caption). The Cδ¹ ¹³C shift ~131 ppm vs Cε¹
~129 ppm vs Cζ ~128 ppm vs Cγ ~138 ppm (Cavanagh, Fairbrother,
Palmer, Skelton 2007 "Protein NMR Spectroscopy" 2e ch. 3) is
recoverable from the typed (Locant, BranchAddress, RingPositionLabel)
tuple.

### Tyr — `RingSystemKind::Benzene_Tyr` (chain `kTyrAtoms`)

**Same layout as Phe + the para-OH atoms.** Atoms `OH`, `HH`
(Locant=Eta, branch={0,0}) carry `PlanarGroupKind::AromaticHydroxyl`
and `PolarHKind::HydroxylOH_Aromatic` for HH, but their
`RingPosition.primary.ring = NotInRing` (correct: the OH is
attached to the ring at Cζ but not in the ring itself).

The TYM variant (`kTyrAtoms_TYM` at line 627-648) replaces
`PlanarGroupKind::AromaticHydroxyl` with `PlanarGroupKind::AromaticOxide`
on OH (deprotonated phenolate, formal_charge=-1 per line 639);
the ring-position labels on the carbons are unchanged. The
phenolate's stronger conjugation with the ring π system (per
`SemanticEnums.h:283-298`, citing Vollhardt & Schore 2018 ch. 22 +
Bovey 1988 ch. 4) is encoded only via the planar-group axis, NOT
via any per-ring label change. This is an open question for
`Distributed ring current` planned-calculator B4 — see Open
Questions §1.

**¹³C shift recovery.** Cδ ~133 (ortho), Cε ~118 (meta — strong
ortho-OH effect), Cζ ~157 (para-OH), Cγ ~127 (ipso) — per
Cavanagh et al. 2007. The (Element, Locant, BranchAddress,
RingPositionLabel) tuple distinguishes these positions; per-position
stratification is reachable.

### His — `RingSystemKind::Imidazole_His` (variants HID/HIE/HIP)

The substrate has THREE separate residue tables — `kHisAtoms_HID`,
`kHisAtoms_HIE`, `kHisAtoms_HIP` (per
`src/generated/LegacyAmberSemanticTables.cpp:446-507`) — selected
at LookupBy time by `Residue.protonation_variant_index`. The
"default HIS" without a variant resolves to HIE per AMBER ff14SB
convention (per `topology-residue-reference-2026-05-05.md` §HIS).

Substrate position labels per atom across the three variants:

| AMBER name | identity | HID position | HIE position | HIP position |
|---|---|---|---|---|
| `CG` | `(C, Gamma, {0,0}, None, None)` | `Ipso/5/t/2` | `Ipso/5/t/2` | `Ipso/5/t/2` |
| `ND1` | `(N, Delta, {1,0}, None, None)` | `Heteroatom_NH/5/t/2` | `Heteroatom_NoH/5/t/2` | `Heteroatom_NH/5/t/2` |
| `CD2` | `(C, Delta, {2,0}, None, None)` | `PyrroleBeta/5/t/2` | `PyrroleBeta/5/t/2` | `PyrroleBeta/5/t/2` |
| `CE1` | `(C, Epsilon, {1,0}, None, None)` | `PyrroleAlpha/5/t/2` | `PyrroleAlpha/5/t/2` | `PyrroleAlpha/5/t/2` |
| `NE2` | `(N, Epsilon, {2,0}, None, None)` | `Heteroatom_NoH/5/t/2` | `Heteroatom_NH/5/t/2` | `Heteroatom_NH/5/t/2` (formal_chg=+1) |

The HID/HIE/HIP per-atom delta is exactly as documented in the
README §5: which N is `Heteroatom_NH` (carries an attached H) vs
`Heteroatom_NoH` (bare N). HIP carries +1 formal charge on Nε2 per
the Lewis-localised convention from
`spec/plan/topology-encoding-dependencies-2026-05-05.md` §C.4.

**Cyclic walk for HIS** in `AminoAcidType.cpp:92`:
`{"CG","ND1","CE1","NE2","CD2"}` — i.e. `Ipso → Heteroatom_X →
PyrroleAlpha → Heteroatom_Y → PyrroleBeta` traversal, where X/Y
flip per variant. This walk is NOT recoverable from
`RingPositionLabel` enum integer order alone: enum walk would be
`Ipso(1) → PyrroleAlpha(7) → PyrroleBeta(8) → Heteroatom_NH(10)
→ Heteroatom_NoH(11)`, which is not the cyclic walk. Bundle C's
`OrderedAtomIndices` for imidazole must encode the chemistry
walk explicitly.

**Tautomer chemistry recovery.** HID Hε1 ~7.4 ppm; HIE Hδ2 ~7.2
ppm; HIP both Hδ2 + Hε1 ~8.5-9.0 ppm (Cavanagh et al. 2007 +
Bertini, Luchinat, Parigi 2001). HID Cδ2 ~120 ppm; HIE Cδ2 ~135
ppm — 15 ppm tautomer-distinguishing signal. Per-tautomer
distinction is reachable from the typed (variant, atom_identity)
join.

### TRP — `RingSystemKind::Indole_Trp_5` and `Indole_Trp_6`

**Per-atom substrate (from `kTrpAtoms`,
`src/generated/LegacyAmberSemanticTables.cpp:373-398`):**

| AMBER name | identity | primary | secondary |
|---|---|---|---|
| `CG` | `(C, Gamma, {0,0}, None, None)` | `Indole_Trp_5 / Ipso / 5 / t / het=1` | `NotInRing` |
| `CD1` | `(C, Delta, {1,0}, None, None)` | `Indole_Trp_5 / PyrroleBeta / 5 / t / 1` | `NotInRing` |
| `NE1` | `(N, Epsilon, {1,0}, None, None)` | `Indole_Trp_5 / Heteroatom_NH / 5 / t / 1` | `NotInRing` |
| `CD2` | `(C, Delta, {2,0}, None, None)` | `Indole_Trp_5 / BridgeFusion / 5 / t / 1` | `Indole_Trp_6 / BridgeFusion / 6 / t / 0` |
| `CE2` | `(C, Epsilon, {2,0}, None, None)` | `Indole_Trp_5 / BridgeFusion / 5 / t / 1` | `Indole_Trp_6 / BridgeFusion / 6 / t / 0` |
| `CE3` | `(C, Epsilon, {3,0}, None, None)` | `Indole_Trp_6 / Ortho1 / 6 / t / 0` | `NotInRing` |
| `CZ2` | `(C, Zeta, {2,0}, None, None)` | `Indole_Trp_6 / Ortho2 / 6 / t / 0` | `NotInRing` |
| `CZ3` | `(C, Zeta, {3,0}, None, None)` | `Indole_Trp_6 / Meta1 / 6 / t / 0` | `NotInRing` |
| `CH2` | `(C, Eta, {2,0}, None, None)` | `Indole_Trp_6 / Meta2 / 6 / t / 0` | `NotInRing` |

The bridgeheads `CD2` and `CE2` are the only atoms with a populated
secondary ring; both are `BridgeFusion` in their primary 5-ring
and `BridgeFusion` (also) in their secondary 6-ring. **The
substrate does not assign a non-bridge label to CD2/CE2 in the
6-ring.** Per the inventory's TRP-6 cyclic walk
(`AminoAcidType.cpp:186`: `{"CD2","CE2","CZ2","CH2","CZ3","CE3"}`),
CD2 is the start atom of the 6-ring walk and CE2 is the next; both
are bridgeheads so the substrate's `BridgeFusion / BridgeFusion`
choice is consistent, but the cyclic-walk ordering is not derivable
from the labels alone.

**TRP-9 perimeter (synthetic ring) — NOT in substrate.** No
`RingSystemKind::Indole_Trp_9` exists; no per-atom secondary
RingPosition labels the perimeter walk. The 9-atom perimeter
ordering today is hand-listed at `AminoAcidType.cpp:188`:
`{"CG","CD1","NE1","CE2","CZ2","CH2","CZ3","CE3","CD2"}` (i.e.
walk the 5-ring from CG to CE2, traverse the 6-ring perimeter from
CE2 to CD2, return to CG). Bundle C's `ConstructRingsFromSubstrate`
must synthesise the perimeter at runtime from the union; the
cyclic ordering is a runtime convention.

**¹³C shift recovery.** Trp pyrrole: Cδ1 (PyrroleBeta) ~124-126
ppm; Cε2 (BridgeFusion 5-ring) ~137-139; Cδ2 (BridgeFusion 5-ring +
6-ring) ~127-128 — per Cavanagh et al. 2007. Trp benzene: Cε3
~118-120 (Ortho1); Cζ2 ~113-115 (Ortho2); Cζ3 ~121-122 (Meta1);
Cη2 ~124-125 (Meta2). The double bridgehead CD2 has the largest
chemical-shift dispersion sensitivity to environment — substrate
captures this via secondary `BridgeFusion` membership.

### PRO — `RingSystemKind::Pyrrolidine_Pro`

**Per-atom substrate (from `kProAtoms`,
`src/generated/LegacyAmberSemanticTables.cpp:322-337`):**

| AMBER name | identity | primary | aromatic | planar |
|---|---|---|---|---|
| `N` | `(N, None, {0,0}, None, Nitrogen)` | `Pyrrolidine_Pro / Heteroatom_NoH / 5 / f / 1` | false | false |
| `CA` | `(C, None, {0,0}, None, AlphaCarbon)` | `Pyrrolidine_Pro / Saturated / 5 / f / 1` | false | false |
| `CB` | `(C, Beta, {0,0}, None, None)` | `Pyrrolidine_Pro / Saturated / 5 / f / 1` | false | false |
| `CG` | `(C, Gamma, {0,0}, None, None)` | `Pyrrolidine_Pro / Saturated / 5 / f / 1` | false | false |
| `CD` | `(C, Delta, {0,0}, None, None)` | `Pyrrolidine_Pro / Saturated / 5 / f / 1` | false | false |

H atoms (`HA`, `HB2`, `HB3`, `HG2`, `HG3`, `HD2`, `HD3`) carry
`RingPosition = NotInRing` even though they are bonded to ring
atoms. This is consistent with the project's "ring atoms = the
heavy atoms of the ring" convention (PHE H atoms at ring positions
are an exception precisely because the aromatic-ring-current
calculators evaluate H positions inside the ring face — not
applicable to Pro saturated chemistry).

**RingPositionLabel coverage.** Pro N is correctly `Heteroatom_NoH`
(secondary amine, no H attached because Pro has no backbone H).
Pro Cα through Cδ all share `Saturated` (value 13). This is
**coarse**: chemistry distinctions among Cα (chiral centre, in ring
+ on backbone), Cβ (sidechain methylene, attached to puckering
pivot), Cγ (puckering pivot — endo/exo mode), Cδ (sidechain
methylene bonded to N) are all collapsed.

**The locked Bundle C extension** adds five typed labels —
`ProRingNitrogen`, `ProRingAlphaCarbon`, `ProRingBeta`,
`ProRingPuckerPivot`, `ProRingDelta` — which restore per-atom
chemistry distinctions that downstream calculators (puckering
classification per Vega & Boyer 1979 Biopolymers 18, 1797-1809;
Schubert, Buynak, Schweitzer-Stenner 2002) read. The extension is
additive to the existing enum; existing `Saturated` value
preserved for any future non-Pro saturated ring chemistry.

**Cyclic walk per locked decision.** `N → Cα → Cβ → Cγ → Cδ`,
matching the chemistry-walk convention from the README. This is
explicitly NOT the residue-walk-order in `kProAtoms` (which lists N,
CA, C, O, CB, CG, CD), nor the chi-angle order from
`AminoAcidType.cpp:156` (`{"N","CA","CB","CG"}` + `{"CA","CB","CG","CD"}`).
The substrate generator currently emits ring atoms in residue-walk
order; Bundle C's `OrderedAtomIndices` for Pro reorders to the
puckering convention.

**No `RingTypeIndex::ProPyrrolidine` today.** Adding it is locked.
No `ProPyrrolidineRing` subclass; locked addition. `Intensity() = 0`
per Joule & Mills 2010 ch. 7 (saturated heterocycles do not carry
circulating π current); this is the load-bearing chemistry fact.

### Ring atoms ON cap-tables (`kCapNtermCharged`, etc.)

For non-PRO N-terminal residues, the cap delta only affects the
backbone N (`PolarHKind`, `PlanarGroupKind`, `formal_charge`). Per
`ApplyCapDelta` at
`src/generated/LegacyAmberSemanticTables.h:107-130`:
**`ring_position` is preserved from chain row, NOT overridden.**
This handles the Pro NTERM correctly: chain N row carries
`Pyrrolidine_Pro/Heteroatom_NoH/5/f/1`; cap delta does not clobber
it. (The fix is documented at lines 119-126 — codex-review
finding from May 2026.)

For non-Pro residues, chain N row carries `NotInRing/NotInRing` for
the backbone N; cap delta preserves it.

### Cap-only atoms (H1/H2/H3 NTERM, OXT/HXT CTERM)

All cap-only atoms have `RingPosition = NotInRing` in the cap
tables (per the `kCapNtermCharged` etc. tables at lines 655-679).
Correct: NTERM H1/H2/H3 are bonded to the backbone N which itself
is in the Pro ring (when residue=Pro), but the H atoms themselves
are not ring atoms. Same for CTERM OXT.


## Gap analysis

Each gap entry: what's missing, why it matters per chemistry /
planned-calculator citation, and the substrate-side recommendation.
Locked decisions are not gaps and not relitigated.

### Gap 1 — `RingTypeIndex::ProPyrrolidine` not present

**What's missing.** `src/Types.h:175-185` declares `RingTypeIndex`
with 8 values + `Count = 8` sentinel (PheBenzene, TyrPhenol,
TrpBenzene, TrpPyrrole, TrpPerimeter, HisImidazole, HidImidazole,
HieImidazole). `ProPyrrolidine` does not exist.

**Why it matters.** Bundle C's `MapRingSystemToRingTypeIndex`
helper (per `chemistry-question-3-ring-topology-2026-05-05.md`
§5.2) requires a `RingTypeIndex` value to map
`RingSystemKind::Pyrrolidine_Pro → RingTypeIndex::ProPyrrolidine`.
`Ring::CreateRing(RingTypeIndex)` likewise needs a case. The
calibration's per-ring-type feature-vector schema (per
`PATTERNS.md` §22) sizes columns by `RingTypeIndex::Count` —
adding the 9th column gives Pro its own (zero-valued) feature
column, preserving the "every ring chemistry in the substrate gets
a feature column" correspondence.

**Recommendation.** Extend `RingTypeIndex` to 9 values with
`ProPyrrolidine = 8` placed AFTER the eight aromatic types.
`Count = 9`. Add an explicit `kAromaticRingTypeCount = 8` named
constant (also in `Types.h`) so the calculator-internal `< 8`
guards become `< kAromaticRingTypeCount` — the per-aromatic-type
arrays' shape is unchanged (locked); the named constant
documents the boundary. Update `RingTypeName(RingTypeIndex)` at
`Types.h:187-199` to add `case ProPyrrolidine: return "PRO5";`
(or `"PRO"`, matching the `Ring::TypeName()` convention).

Citation: Joule, J. A. & Mills, K. (2010) "Heterocyclic Chemistry,"
5th ed., Wiley, ch. 7 — saturated heterocycles do not carry
circulating π current. The `Intensity()` returning 0 is the
chemistry-load-bearing fact; the index value extension is the
typed-substrate carrier.

### Gap 2 — `ProPyrrolidineRing` class absent from `Ring.h`

**What's missing.** `Ring.h:50-202` declares `Ring`,
`SixMemberedRing`, `FiveMemberedRing`, `FusedRing`, plus eight
concrete subclasses (`PheBenzeneRing`, `TyrPhenolRing`,
`TrpBenzeneRing`, `TrpPyrroleRing`, `HisImidazoleRing`,
`HidImidazoleRing`, `HieImidazoleRing`, `IndolePerimeterRing`). No
class for Pro pyrrolidine.

**Why it matters.** The Ring class hierarchy is the typed surface
calculators consume; adding Pro chemistry means adding a class
that returns the Pro chemistry from `Aromaticity()`,
`RingSizeValue()`, `NitrogenCount()`, `Intensity()`, `JBLobeOffset()`,
`LiteratureIntensity()`, `TypeName()`. The Ring object substrate
discipline says objects answer questions about themselves — Pro
pyrrolidine ring's chemistry is encoded in those override methods.

**Recommendation.** Add `ProPyrrolidineRing` as a subclass of
`FiveMemberedRing` (or introduce a `SaturatedRing` base sister to
`SixMemberedRing` / `FiveMemberedRing` / `FusedRing` if the audit
favours consistency with Ring.h's category-by-shape pattern).
Override values:

```text
type_index           = RingTypeIndex::ProPyrrolidine
Intensity()          = CalculatorConfig::Get("pro_pyrrolidine_ring_current_intensity") (or 0.0 directly)
LiteratureIntensity() = 0.0
JBLobeOffset()       = CalculatorConfig::Get("pro_pyrrolidine_jb_lobe_offset") (any reasonable saturated-ring value; not used while Intensity=0)
NitrogenCount()      = 1
Aromaticity()        = RingAromaticity::None  // requires extending RingAromaticity enum (Gap 2a)
RingSizeValue()      = 5
TypeName()           = "PRO5"
```

Cite Joule & Mills 2010 ch. 7 in the docstring (`Intensity()=0`
because saturated heterocycles do not support a circulating π
current).

**Sub-gap 2a — `RingAromaticity::None` value.** `Types.h:168`
declares `enum class RingAromaticity { Full, Reduced, Weak };`.
Saturated rings (Pro pyrrolidine; future heterocycle subclasses)
need a fourth value `None`. Citation: Vollhardt & Schore 2018
ch. 14 (aromaticity definition: cyclic, planar, fully conjugated, 4n+2
π electrons — saturated rings fail all four criteria).

Add `CreateRing(RingTypeIndex::ProPyrrolidine)` case to
`src/Ring.cpp:50-62`:

```cpp
case RingTypeIndex::ProPyrrolidine: return std::make_unique<ProPyrrolidineRing>();
```

The factory's `default:` clause currently returns
`PheBenzeneRing` — that fallback is a soft-defect (a misclassified
or future RingTypeIndex silently becomes a Phe ring). Bundle C is
a natural place to **fail-loudly** instead (per the
substrate-miss FATAL discipline from Bundle B); not strictly
substrate-side, flagged here as adjacent.

### Gap 3 — Five Pro-specific `RingPositionLabel` values absent

**What's missing.** `RingPositionLabel` (`SemanticEnums.h:528-589`)
has 14 values; the only saturated-ring value is `Saturated = 13`.
The Pro N is `Heteroatom_NoH` (value 11 — generic to bare N in
any ring); Pro Cα/Cβ/Cγ/Cδ all collapse to `Saturated`. Per the
locked decision, five typed Pro-specific labels are added.

**Why it matters.** The chemistry distinctions among Pro ring
atoms are real and load-bearing:

- **Cα** is on the backbone (BackboneRole=AlphaCarbon) AND in the
  ring (RingSystemKind=Pyrrolidine_Pro). It's a chiral centre
  participating in puckering geometry. (Vega & Boyer 1979.)
- **Cβ** is the side-chain methylene attached to the chiral Cα.
  Cβ ¹³C exo-vs-endo shift difference is ~1-2 ppm (Sarkar, Young,
  Sullivan 1986; Schubert, Buynak, Schweitzer-Stenner 2002).
- **Cγ** is the puckering pivot — the atom that moves above/below
  the Cα-Cδ-N plane in endo vs exo. Its ¹³C shift difference between
  puckers is the largest of the five (Schubert et al. 2002).
- **Cδ** is the side-chain methylene bonded to the ring N. The
  Cδ-N bond is what closes the ring and what dictates the
  trans/cis isomerism of the preceding peptide bond.
- **N** is the secondary amine carrying the (i-1) C-N peptide
  bond. Pro's lack of backbone H is precisely because N is in the
  ring as a secondary amine.

Future calculators reading `ring_position.primary.position ==
ProRingPuckerPivot` recover Cγ chemistry directly without
recovering it from `(Locant=Gamma + BackboneRole=None +
ring_membership=Pyrrolidine_Pro)` — the typed surface IS the
chemistry.

**Recommendation.** Per the locked decision, extend
`RingPositionLabel` (additive, append to enum after `Saturated`):

```text
ProRingNitrogen       = 14   ///< Pro N (backbone, secondary amine, in ring; absent backbone H).
ProRingAlphaCarbon    = 15   ///< Pro Cα (backbone, chiral centre, in ring).
ProRingBeta           = 16   ///< Pro Cβ (sidechain methylene, attached to chiral Cα).
ProRingPuckerPivot    = 17   ///< Pro Cγ (puckering pivot — endo/exo mode).
ProRingDelta          = 18   ///< Pro Cδ (sidechain methylene bonded to N).
```

Existing `Saturated = 13` value is preserved (Bundle C does not
remove enum values). Update Pro substrate-table generator
output (`tools/topology/build_semantic_tables.cpp`) to emit the
new labels for Pro-residue ring atoms; regenerate
`src/generated/LegacyAmberSemanticTables.cpp`. Existing 30 residue
tables (20 standard + 10 variants) plus 4 cap tables continue to
emit `Saturated` for any non-Pro saturated ring (none exist in the
standard 20, but the value remains for future extensions).

Update reference doc
`spec/plan/topology-residue-reference-2026-05-05.md` PRO section
(lines 609-634) to use the new labels. The existing
`Pyrrolidine_Pro/Saturated/5/f/1` rows for Cα through Cδ become
`Pyrrolidine_Pro/ProRingAlphaCarbon/5/f/1`,
`Pyrrolidine_Pro/ProRingBeta/5/f/1`, etc. The existing
`Pyrrolidine_Pro/Heteroatom_NoH/5/f/1` row for N becomes
`Pyrrolidine_Pro/ProRingNitrogen/5/f/1`.

Citations: Vega & Boyer (1979) "Conformational analysis of proline
residues in a peptide chain," Biopolymers 18, 1797-1809 (puckering
chemistry); Sarkar, Young, Sullivan (1986) (¹³C shift dependence);
Schubert, Buynak, Schweitzer-Stenner (2002) (modern classification
methodology); Joule & Mills 2010 ch. 7 (saturated 5-ring chemistry).

### Gap 4 — TRP perimeter has no substrate encoding

**What's missing.** The 9-atom TRP indole perimeter
(`RingTypeIndex::TrpPerimeter`) has no corresponding
`RingSystemKind` value; no per-atom secondary `RingPosition`
labels mark perimeter ortho/meta/para. The 9-atom union and
ordering are runtime synthesis constructs in
`Protein::DetectAromaticRings` today (and in
`AminoAcidType.cpp:188`).

**Why it matters.** Per the README §6 (citing Case 1995, J. Biomol.
NMR 6, 341-346 + PATTERNS.md §24): the conjugated π current in
indole flows around all 9 atoms as a single circuit — modeling
TRP as just `TrpPyrrole` + `TrpBenzene` undercounts the perimeter
contribution. The empirical I = -19.2 nA/T = I(TrpPyrrole) +
I(TrpBenzene) is verified at 1.000 ratio post-2026-04-02 (HM
correction). The perimeter ring is a calculation construct, but
calculators consume it as a Ring object today and must continue
to. Bundle C must materialise it.

**Two options for substrate-side encoding:**

- **Option (a) — Synthesise at runtime in `ConstructRingsFromSubstrate`.**
  No substrate change; perimeter atom-list = union of `Indole_Trp_5`
  atoms + `Indole_Trp_6` atoms minus shared bridgeheads, with a
  defined cyclic order. The substrate continues to encode only
  the chemical rings (`Indole_Trp_5` and `Indole_Trp_6`).

- **Option (b) — Extend substrate with `RingSystemKind::Indole_Trp_9`.**
  Add the 9th `RingSystemKind` value; substrate generator emits
  per-atom secondary RingMembership records for TRP perimeter
  atoms; the perimeter cyclic walk encoded via additional
  `RingPositionLabel` values like `PerimeterOrtho1/Ortho2/...`.
  Per-position-on-perimeter calculators read perimeter labels
  directly.

**Recommendation.** Per the README's note on TRP perimeter
representation as "substrate-extension as strongly recommended
default" — **Option (b) is the recommended substrate posture**.
Practical considerations:

- Adding `RingSystemKind::Indole_Trp_9` to `SemanticEnums.h:496-504`.
- Adding perimeter-walk position labels as additive
  `RingPositionLabel` extensions: e.g. `PerimeterIpso`,
  `PerimeterOrtho1`, `PerimeterOrtho2`, `PerimeterMeta1`,
  `PerimeterMeta2`, `PerimeterPara`, `PerimeterBridgeIn`,
  `PerimeterBridgeOut`, `PerimeterPyrroleNitrogen` (one per
  perimeter atom; or coarser, share with existing labels by axis).
  This is more enum work; alternative is to use existing labels
  (`Ipso → BridgeFusion → Heteroatom_NH → BridgeFusion → Ortho1
  → Ortho2 → Meta1 → Meta2 → BridgeFusion`) at the cost of
  reusing labels in non-natural ways.
- Substrate generator emits TRP perimeter atoms with secondary
  `RingPosition` populated. The 5-ring's primary positions and
  the 6-ring's primary positions are unchanged; only the secondary
  field on the 9 atoms gains perimeter membership.

The audit notes that Option (a) is operationally simpler and the
chemistry's stable enough that runtime synthesis is not a
maintenance burden. The user's locked-decision note ("substrate
extension as strongly recommended default; permission to argue
runtime synthesis if a per-calculator agent has a science case")
suggests Option (b) unless a calculator has an objection.
**The audit's recommendation: Option (b)**, with the caveat that
the new `RingPositionLabel` values be chosen carefully to avoid
overlap with the existing per-atom labels of `Indole_Trp_5` and
`Indole_Trp_6` (each Trp atom would carry primary in its chemical
ring + secondary in the perimeter ring + the perimeter ring would
have its own walk-position labels distinct from primary).

If Bundle C ships Option (a) for sequencing reasons (less
substrate-generator work, scoped tightly), the gap remains for
Phase E or a later substrate iteration; the cyclic-ordering
convention for the perimeter walk needs to be documented in
`Ring.h` regardless.

Citations: Case (1995) "Calibration of ring-current models for the
calculation of protein chemical shifts," J. Biomol. NMR 6, 341-346
(perimeter intensity calibration); Vollhardt & Schore 2018 ch. 16
(fused-aromatic-ring conjugation chemistry).

### Gap 5 — `Ring::atom_indices` ordering convention undocumented

**What's missing.** The cyclic walk through ring atoms determines
the SVD ring-normal sign (per `Ring.cpp:30-37`; verified
analytically per PATTERNS.md §22 — Phe at I = -12.0, probe 3 Å
axial → σ = +1.40 ppm shielded, sign convention catches a bug
not caught by compilation or unit tests). The convention is
load-bearing for ring-current calculator output; the substrate
must drive a deterministic walk.

Today's de facto convention (from `AminoAcidType.cpp` ring atom
lists):

```text
PHE: CG → CD1 → CE1 → CZ → CE2 → CD2
TYR: CG → CD1 → CE1 → CZ → CE2 → CD2
TRP6: CD2 → CE2 → CZ2 → CH2 → CZ3 → CE3
TRP5: CG → CD1 → NE1 → CE2 → CD2
TRP9: CG → CD1 → NE1 → CE2 → CZ2 → CH2 → CZ3 → CE3 → CD2
HIS: CG → ND1 → CE1 → NE2 → CD2
```

In typed RingPositionLabel terms:

```text
PHE/TYR: Ipso → Ortho1 → Meta1 → Para → Meta2 → Ortho2  (6-ring benzene walk)
TRP6: BridgeFusion(start) → BridgeFusion → Ortho2 → Meta2 → Meta1 → Ortho1
TRP5: Ipso → PyrroleBeta → Heteroatom_NH → BridgeFusion → BridgeFusion
HIS: Ipso → Heteroatom_X → PyrroleAlpha → Heteroatom_Y → PyrroleBeta
```

Notes per chemistry:

- The PHE / TYR walk Ipso → Ortho1 → Meta1 → Para → Meta2 → Ortho2
  is the natural cyclic walk around the 6-ring with Ortho1/Ortho2
  branched away from the ipso → para line. The substrate's
  Ortho1/Ortho2 distinguishes CD1/CD2 by the Markley 1998
  chi₂-priority convention (`SemanticEnums.h:540`).
- TRP-6's start-atom `CD2` is a BridgeFusion (shared with TRP-5),
  not an Ipso. Substrate doesn't surface TRP-6's "ipso"
  conceptually — the shared bridge IS the boundary; the walk
  starts at one bridge and ends at the other.
- The Pro convention `N → Cα → Cβ → Cγ → Cδ` is the locked
  Bundle C choice — the residue-walk-order ordering, NOT the
  geometric cyclic order (the cyclic order is `N → Cδ → Cγ → Cβ
  → Cα` going one way around, `N → Cα → Cβ → Cγ → Cδ` going the
  other way; both close the 5-ring). The locked choice is
  `N → Cα → Cβ → Cγ → Cδ` matching the residue's natural
  outward-from-backbone walk.

**Why it matters.** Bundle C's `OrderedAtomIndices(atom_list,
RingSystemKind, Residue)` helper walks atoms in canonical typed
order. The convention chosen MUST match today's de facto order
bit-identically for non-Pro residues, OR a deliberate sign-flip
must be documented + the SVD-normal sign convention re-verified
on bless fixtures. Per `chemistry-question-3-ring-topology-2026-05-05.md`
§7.3, this is a verified-then-locked decision.

**Recommendation.** Document the cyclic-walk convention per
`RingTypeIndex` in `Ring.h` as a docstring above each subclass
declaration. For example:

```cpp
class PheBenzeneRing : public SixMemberedRing {
public:
    // atom_indices cyclic walk: Ipso → Ortho1 → Meta1 → Para → Meta2 → Ortho2.
    // Heavy atoms in Markley 1998 §2.1.1 convention (chi₂-priority for
    // Ortho1 vs Ortho2). Ring-normal sign convention: cross product of
    // first two edges (Ring::ComputeGeometry, src/Ring.cpp:34-37).
    PheBenzeneRing() { type_index = RingTypeIndex::PheBenzene; }
    ...
};
```

Repeat for the eight subclasses + ProPyrrolidineRing. Bundle C's
`OrderedAtomIndices` reads these conventions out of the typed
labels (the docstring documents what code computes). Add a unit
test that compares post-Bundle-C cyclic ordering against the
pre-Bundle-C ordering on a canonical fixture (1UBQ + 1Z9B + 1P9J)
and asserts (residue, type_index, atom_indices in cyclic order)
match exactly. Per
`chemistry-question-3-ring-topology-2026-05-05.md` §5.6.

### Gap 6 — HIP imidazolium charge asymmetry not surfaced

**What's missing.** HIP carries +1 formal charge, distributed
asymmetrically per the Lewis-localised convention: per
`spec/plan/topology-encoding-dependencies-2026-05-05.md` §C.4, the
+1 sits on Nε2 (per the kHisAtoms_HIP table line 498:
`formal_charge = 1` on NE2), Nδ1 has formal_charge = 0. Both
nitrogens are tagged `Heteroatom_NH` in the substrate
(`kHisAtoms_HIP` entries 6 + 9, plus the H atoms HD1 + HE2 carry
`PolarHKind::ImidazoleNH` per lines 503 + 506).

**Why it matters.** The imidazolium charge asymmetry, even though
delocalised at the electron-density level, IS chemistry — the
Lewis-localised typed encoding distinguishes which nitrogen bears
the formal charge. Per the README §5:

```text
HIP both Hδ2 + Hε1 ~8.5-9.0 ppm (charged ring shifts ring protons strongly downfield)
```

Calculators that gate on formal_charge (planned-audit F13:
A11 charged-residue strata) need to find the charge on Nε2
specifically.

**Status: not a gap requiring substrate extension.** The Lewis
formal charge on Nε2 is correctly encoded at `kHisAtoms_HIP` line
498 (`formal_charge = 1`); the chemistry IS represented. The
finding here is documentary — calculator authors reading
`SemanticAt(NE2_atom_idx).formal_charge` find +1 (HIP only); on
HID/HIE the same lookup returns 0. The asymmetry is reachable.

The audit flags this only because the per-atom RingPositionLabel
on both nitrogens is `Heteroatom_NH` for HIP, which loses the
"this N is the charged one" axis at the ring-position-label
granularity. A future fine-grained label distinction (e.g.
`Heteroatom_NH_Charged` vs `Heteroatom_NH_Neutral`) is possible
but may be redundant with `formal_charge`.

**Recommendation.** No change — the chemistry is reachable via
(`RingPositionLabel == Heteroatom_NH` AND `formal_charge != 0`)
join. Document this access path in `Ring.h` or the ring-construction
docstring.

### Gap 7 — No residue-level "which rings" predicate

**What's missing.** `LegacyAmberTopology` does not expose
`std::vector<RingSystemKind> RingsForResidue(size_t residue_index)`.
Today's runtime code computes this by iterating per-atom
`SemanticAt(ai).ring_position.primary.ring` and collecting unique
values per residue. The current `Protein::DetectAromaticRings`
iterates `aatype.rings[].atom_names` instead — a string table.
Bundle C's `ConstructRingsFromSubstrate` will iterate per-atom
substrate; the residue-level groupings happen there.

**Why it matters.** This is convenience: a UI / debug consumer
that wants "what rings does residue 47 have?" without going
through Protein::RingAt iteration finds no typed predicate
today. Planned-calculator audit's TS8 (Aromatic-H sanity check)
joins on this answer. The substrate has the data; the
interrogation surface is by-atom not by-residue.

**Recommendation.** Optional: add to `LegacyAmberTopology`:

```cpp
// Returns the unique RingSystemKinds that any atom in this residue
// participates in (primary or secondary). Empty vector for
// non-ring residues.
std::vector<RingSystemKind> RingsForResidue(
    size_t residue_index,
    const std::vector<Residue>& residues) const;
```

Implementation iterates the residue's atom_indices, collects unique
`SemanticAt(ai).ring_position.primary.ring` values (and secondary
where non-NotInRing). This is utility; not essential for Bundle C
(which iterates by-atom directly per the §5.2 design).

### Gap 8 — No typed RingPositionLabel → IUPAC / BMRB / AMBER atom name projection

**What's missing.** Session E (per the project's roadmap) lands
projection functions on `LegacyAmberTopology`: `IupacNameAt`,
`BmrbNameAt`, `AmberNameAt` from typed substrate. For ring atoms,
these projections are nearly identity functions:

- AMBER name lives directly on the runtime `Atom.pdb_atom_name`
  (today's storage; canonical post-NamingApplicator).
- IUPAC name for ring atoms = AMBER name (modulo Markley 1998
  `δ` Greek-letter rendering vs AMBER `D` Roman rendering — same
  storage integer).
- BMRB nomenclature page at `https://bmrb.io/referenc/nomenclature`
  follows IUPAC for ring atoms in the standard 20.

For HID/HIE/HIP the projection has to read
`Residue.protonation_variant_index` to decide which variant the
atom is from; the typed substrate already does that via LookupBy.

**Why it matters.** Projection-readiness for ring atoms is shallow
because BMRB ≡ IUPAC ≡ AMBER for ring atoms in the standard 20.
The audit confirms substrate is sufficient: every ring atom is
reachable via `(residue index, AtomMechanicalIdentity)` lookup
which selects the typed substrate row, which has the atom's
identity (Element, Locant, BranchAddress, etc.). Re-deriving the
name from the typed identity is straightforward.

**Recommendation.** Session E (after Bundle C) lands the
projection functions. The substrate is ready for them as-is —
no Bundle C-time substrate extension needed for projection. Flag
for Session E sequencing: the `RingPositionLabel` itself is NOT
the atom name; it is the chemistry-position label. The atom name
comes from `Locant + BranchAddress + DiastereotopicIndex +
BackboneRole` (the existing AtomMechanicalIdentity surface).
Session E's projection functions use that.

For Pro post-Bundle-C, the new `ProRingPuckerPivot` etc. labels
do NOT enter the IUPAC / BMRB / AMBER projection — Cγ's atom name
is "CG" via Locant=Gamma, not via the new ring-position label.
Projection reads the identity tuple. Confirms: gap 8 stays "no
substrate work needed for projection."

### Gap 9 — `Ring` class lacks per-vertex RingPositionLabel exposure

**What's missing.** `Ring::atom_indices` is `std::vector<size_t>`;
calculators iterating ring atoms get atom indices but no typed
position label. To recover position for vertex k, calculator must
call `protein.LegacyAmber().SemanticAt(ring.atom_indices[k]).ring_position.primary.position`
— the access path works (per planned-audit §B Constraint B), but
the Ring object itself does not expose a per-vertex label
accessor.

**Why it matters.** Phase 2 / planned-calculator §B4 (Distributed
ring current per planned-audit lines 264-272) needs per-vertex
labels to decompose contributions per ring position. The audit
expression "constraint B: substrate already exposes ring_position
per atom via SemanticAt(ai). The Ring object itself does NOT need
a per-vertex `RingPositionLabel` member (redundant with substrate).
Calculators that want per-position info read substrate via
atom_indices." — this is a calculator-side reading discipline,
already accessible.

**Status: not a gap.** The substrate IS reachable via the existing
access path. The Ring class does not need extension; the per-vertex
label is one substrate query away.

The audit notes this only to be explicit: a future bundle that
chooses to add `Ring::PositionLabelAt(k)` would be convenience,
not chemistry. The architectural rule (per
`chemistry-question-3-ring-topology-2026-05-05.md` §6.2) prefers
calculators to read by-atom from substrate when the question is
per-position; Ring class stays minimal.


## Recommended substrate extensions

Concrete code surfaces that need additive change to land Bundle C
in its locked shape. Substrate-side scope only.

### 1. Type enum extensions

**`src/Types.h:175-185`** — extend `RingTypeIndex`:

```text
enum class RingTypeIndex {
    PheBenzene    = 0,
    TyrPhenol     = 1,
    TrpBenzene    = 2,
    TrpPyrrole    = 3,
    TrpPerimeter  = 4,
    HisImidazole  = 5,
    HidImidazole  = 6,
    HieImidazole  = 7,
    ProPyrrolidine = 8,            // ADDED — Bundle C
    Count         = 9              // bumped from 8 to 9
};
inline constexpr int kAromaticRingTypeCount = 8;   // ADDED — names the per-aromatic boundary
```

Update `RingTypeName(RingTypeIndex)` switch to add
`case RingTypeIndex::ProPyrrolidine: return "PRO5";` (or `"PRO"`).

**`src/Types.h:168`** — extend `RingAromaticity`:

```text
enum class RingAromaticity { Full, Reduced, Weak, None };
```

Cite Vollhardt & Schore 2018 ch. 14 in the docstring (saturated
rings fail aromaticity criteria: cyclic + planar + fully conjugated +
4n+2 π electrons).

### 2. SemanticEnums extensions

**`src/SemanticEnums.h:528-589`** — extend `RingPositionLabel`
(append, preserve existing values):

```text
enum class RingPositionLabel : uint8_t {
    NotInRing       = 0,
    Ipso            = 1,
    Ortho1          = 2,
    Ortho2          = 3,
    Meta1           = 4,
    Meta2           = 5,
    Para            = 6,
    PyrroleAlpha    = 7,
    PyrroleBeta     = 8,
    BridgeFusion    = 9,
    Heteroatom_NH   = 10,
    Heteroatom_NoH  = 11,
    Heteroatom_OH   = 12,
    Saturated       = 13,         // PRESERVED — future non-Pro saturated rings
    ProRingNitrogen     = 14,     // ADDED — Pro N (sec amine; in ring; no H)
    ProRingAlphaCarbon  = 15,     // ADDED — Pro Cα (backbone, chiral, in ring)
    ProRingBeta         = 16,     // ADDED — Pro Cβ (sidechain methylene)
    ProRingPuckerPivot  = 17,     // ADDED — Pro Cγ (puckering pivot — endo/exo)
    ProRingDelta        = 18,     // ADDED — Pro Cδ (sidechain methylene bonded to N)
};
```

Each value's docstring cites: Vega & Boyer (1979) Biopolymers 18,
1797-1809; Schubert, Buynak, Schweitzer-Stenner (2002); Joule &
Mills 2010 ch. 7.

**Optional — TRP-9 perimeter substrate encoding (§Gap 4 Option (b)):**

If Option (b) chosen, additive RingSystemKind:

```text
enum class RingSystemKind : uint8_t {
    NotInRing       = 0,
    Benzene_Phe     = 1,
    Benzene_Tyr     = 2,
    Imidazole_His   = 3,
    Indole_Trp_5    = 4,
    Indole_Trp_6    = 5,
    Pyrrolidine_Pro = 6,
    Indole_Trp_9    = 7,           // ADDED — perimeter (synthetic, Case 1995)
};
```

And RingPositionLabel additions for perimeter walk (audit
recommends specific labels per the cyclic-walk convention; flagged
in Open Questions §2 for user blessing).

### 3. Ring class hierarchy extension

**`src/Ring.h`** — add `ProPyrrolidineRing`:

```cpp
class ProPyrrolidineRing : public FiveMemberedRing {
public:
    ProPyrrolidineRing() { type_index = RingTypeIndex::ProPyrrolidine; }
    // No ring current: saturated heterocycle (Joule & Mills 2010 ch. 7).
    double Intensity() const override { return CalculatorConfig::Get("pro_pyrrolidine_ring_current_intensity"); }
    double LiteratureIntensity() const override { return 0.0; }
    // Lobe offset is unused while Intensity = 0; provided for completeness
    // and any future puckering-aware calculator that wants the geometry constant.
    double JBLobeOffset() const override { return CalculatorConfig::Get("pro_pyrrolidine_jb_lobe_offset"); }
    int NitrogenCount() const override { return 1; }
    RingAromaticity Aromaticity() const override { return RingAromaticity::None; }
    // RingSizeValue inherited from FiveMemberedRing (returns 5).
    const char* TypeName() const override { return "PRO5"; }
    // Cyclic walk: N → Cα → Cβ → Cγ → Cδ (locked Bundle C convention,
    // matching residue-walk-order outward from the backbone N).
};
```

**`src/Ring.cpp:50-62`** — extend `CreateRing(RingTypeIndex)`:

```cpp
case RingTypeIndex::ProPyrrolidine: return std::make_unique<ProPyrrolidineRing>();
```

The `default:` fallback to PheBenzene at line 60-61 should be
changed to `assert(false && "Unknown RingTypeIndex"); abort();` per
the substrate-miss FATAL discipline; flagged as adjacent (not
substrate-side).

### 4. Substrate generator regeneration

**`tools/topology/build_semantic_tables.cpp`** — emit the five new
`RingPositionLabel` values at Pro-residue ring atoms; regenerate
`src/generated/LegacyAmberSemanticTables.cpp`. Pre/post comparison:

Before (Pro N, line 323):
```text
{nmr::RingSystemKind::Pyrrolidine_Pro, nmr::RingPositionLabel::Heteroatom_NoH, 5, false, false, 1}
```

After:
```text
{nmr::RingSystemKind::Pyrrolidine_Pro, nmr::RingPositionLabel::ProRingNitrogen, 5, false, false, 1}
```

Same for Cα, Cβ, Cγ, Cδ (changing from `Saturated` to the
respective new labels).

If Option (b) for TRP perimeter chosen, generator also emits
secondary RingPosition records on TRP atoms (CG, CD1, NE1, CE2,
CZ2, CH2, CZ3, CE3, CD2 — 9 atoms total) carrying `Indole_Trp_9`
membership and the 9-atom-walk position labels.

### 5. Reference doc updates

**`spec/plan/topology-residue-reference-2026-05-05.md` §PRO
(lines 609-634)** — update the ring_primary column for Pro
ring atoms:

```text
| N | N | None | NA | — | NotPolar | Pyrrolidine_Pro/ProRingNitrogen/5/f/1 | NotInRing | NotProchiral | 0 | f | secondary amine; in the ring; no H |
| CA | C | None | NA | — | NotPolar | Pyrrolidine_Pro/ProRingAlphaCarbon/5/f/1 | NotInRing | NotProchiral | 0 | f | Cα also in ring; chiral; pucker geometry |
| CB | C | None | NA | — | NotPolar | Pyrrolidine_Pro/ProRingBeta/5/f/1 | NotInRing | NotProchiral | 0 | f | sidechain methylene |
| CG | C | None | NA | — | NotPolar | Pyrrolidine_Pro/ProRingPuckerPivot/5/f/1 | NotInRing | NotProchiral | 0 | f | puckering pivot — endo/exo (Schubert 2002) |
| CD | C | None | NA | — | NotPolar | Pyrrolidine_Pro/ProRingDelta/5/f/1 | NotInRing | NotProchiral | 0 | f | sidechain methylene; bonded to N |
```

**`spec/plan/topology-encoding-dependencies-2026-05-05.md`** — if
this document references `RingPositionLabel::Saturated` for Pro
specifically (audit didn't grep this; flag for documentation
review during Bundle C drafting).

### 6. Ring docstring updates

**`src/Ring.h`** — add cyclic-walk convention docstring above each
ring subclass declaration (per Gap 5). Document:

- PHE/TYR: Ipso → Ortho1 → Meta1 → Para → Meta2 → Ortho2.
- TRP-6: BridgeFusion(CD2-start) → BridgeFusion(CE2) → Ortho2 →
  Meta2 → Meta1 → Ortho1 (residue-walk-order from `AminoAcidType.cpp`).
- TRP-5: Ipso(CG) → PyrroleBeta(CD1) → Heteroatom_NH(NE1) →
  BridgeFusion(CE2) → BridgeFusion(CD2).
- TRP-9 perimeter: CG → CD1 → NE1 → CE2 → CZ2 → CH2 → CZ3 → CE3 → CD2
  (synthetic walk; runtime-defined or substrate-emitted per Option
  choice).
- HID/HIE/HIP: Ipso(CG) → Heteroatom_X(ND1) → PyrroleAlpha(CE1) →
  Heteroatom_Y(NE2) → PyrroleBeta(CD2). X/Y per variant: HID →
  X=NH, Y=NoH; HIE → X=NoH, Y=NH; HIP → X=NH, Y=NH(charged).
- Pro: N → Cα → Cβ → Cγ → Cδ (locked).

Each convention's docstring cites the literature source for the
atom labelling (Markley 1998 §2.1.1 for the side-chain locants;
Joule & Mills 2010 ch. 13 for pyrrole positions; Vega & Boyer 1979
for Pro).


## Categorical topology recommendations

These are typed structural-relationship surfaces the substrate
could expose to enable richer downstream science without
introducing strings. None is required for Bundle C in its locked
shape; all are flagged for future substrate iterations.

### 1. Per-residue ring catalogue accessor

Per Gap 7 — `LegacyAmberTopology::RingsForResidue(residue_index,
residues)` returning the unique `RingSystemKind`s any atom in the
residue participates in (primary or secondary). Implementation
trivial; documents the typed access path as the canonical query.

### 2. Per-residue ring atom set accessor

`LegacyAmberTopology::RingAtomsInResidue(residue_index, ring_kind,
residues)` returning the atom indices (in the residue) whose
primary or secondary `ring_position.ring == ring_kind`. Bundle C's
`ConstructRingsFromSubstrate` does this internally; exposing it
makes per-residue ring-atom enumeration available to UI / debug
without re-iterating.

### 3. Ring-atom adjacency relation (within ring)

Today, calculators that need ring-atom neighbours
(within-ring edges) walk `Ring::atom_indices[k]` and `[k+1 mod N]`,
relying on the cyclic-walk convention. A typed predicate
`AreRingNeighbours(atom_a, atom_b)` reading from substrate (both
atoms have non-NotInRing primary, primary.ring matches, neighbours
in the cyclic order) makes the relation queryable directly. Lower
priority — calculator-internal iteration is the established
pattern.

### 4. Fused-ring relation surface

`Ring::IsFused()` exists at `Ring.h:71`, returning true when
`fused_partner_index != SIZE_MAX`. The substrate-side complement
is per-atom: `ring_position.secondary.ring != NotInRing` indicates
the atom is at a fused-ring bridge. Today this is reachable. A
typed `LegacyAmberTopology::IsFusedRingBridge(atom_index)` accessor
would be utility; not essential.

### 5. Heteroatom-position-class predicate

A predicate that classifies a ring atom's heteroatom relationship:
ring carbon, ring nitrogen-with-H, ring nitrogen-without-H, ring
oxygen, etc. — already encoded in `RingPositionLabel`'s
`Heteroatom_NH` / `Heteroatom_NoH` / `Heteroatom_OH` values.
Calculators that want "all aromatic carbons" gate on
`ring_position.primary.aromatic == true && element == C`. Today
this is reachable; the predicate is utility.

### 6. Ring-membership-only predicate (heavy atom vs H atom)

A predicate that distinguishes ring heavy atoms (which DO carry
ring-current geometry — vertex SVD, etc.) from ring-attached H
atoms (which carry the ring chemistry but not the ring vertex
participation). Today PHE has `H` atoms with non-NotInRing
RingPosition (matching their bonded heavy atom's position) — the
H atoms are tagged as ring-position-bearing for stratified
shielding calculation, not for vertex SVD. The substrate has the
information; the calculator-side discipline reads
`Element != H` to filter.

### 7. n_heteroatoms encoded but underused

`RingMembership.n_heteroatoms` (`SemanticEnums.h:606`) carries the
ring's heteroatom count (0 for Phe/Tyr/Trp-6; 1 for Trp-5/Pro;
2 for HIS imidazole). Planned-calculator audit §B4 mentions this
field for HID/HIE asymmetric per-vertex models; today it's
populated correctly across all generated tables but no consumer
yet uses it. Documented as substrate-ready for future per-ring
heteroatom-count-aware kernels.


## Projection-readiness assessment (BMRB / IUPAC / AMBER per ring atom)

For Session E's projection functions on `LegacyAmberTopology`,
the substrate readiness for ring atoms specifically:

### AMBER atom name

Today on `Atom.pdb_atom_name` (canonical post-NamingApplicator).
Reachable directly. **Substrate-ready.**

For Bundle C purposes, the ring construction path does NOT depend
on `Atom.pdb_atom_name` — the typed substrate
(`AtomMechanicalIdentity` for lookup; `RingPosition` for grouping)
is the chemistry authority. AMBER name is only consumed when
projection is asked for explicitly, which Session E handles.

### IUPAC atom name (per Markley 1998 §2.1.1)

For ring atoms, IUPAC name = AMBER name modulo Greek-letter
rendering (δ vs D, ε vs E, ζ vs Z, η vs H). Storage integer is
the same — the difference is rendering.

The typed projection: `IupacNameAt(atom_index)` reads
`SemanticAt(atom_index).{element, locant, branch, di_index,
backbone_role}` and renders per Markley convention (Greek
letters; italic vs Roman per Markley's text rules; Position2 vs
Position3 disambiguator for diastereotopic methylene Hs).

For Pro post-Bundle-C: same access path. The new
`ProRingPuckerPivot` etc. labels are NOT used in projection; Cγ's
IUPAC name is "Cγ" via Locant=Gamma + Element=C.

**Substrate-ready.**

### BMRB atom name

Per the BMRB nomenclature page
(`https://bmrb.io/referenc/nomenclature`) and pseudoatom table
(`https://bmrb.io/ref_info/pseudoatom_nom.txt`), BMRB names for
ring atoms in the standard 20 align with IUPAC. For HID/HIE/HIP,
BMRB tracks the tautomer-specific Hδ1/Hε2 naming that AMBER
preserves (HID has Hδ1; HIE has Hε2; HIP has both).

The typed projection: `BmrbNameAt(atom_index)` reads the same
(element, locant, branch, di_index, backbone_role) tuple plus
`Residue.protonation_variant_index` for HIS variants. Identical
to IUPAC modulo BMRB's slight rendering preferences (Roman vs
Greek; BMRB uses Roman in deposit format).

**Substrate-ready.**

### Per-position chemistry label (if requested as a name)

The substrate's `RingPositionLabel` is a chemistry-position label,
not an atom name. "Para C of Tyr" maps to atom name "CZ" via
typed reasoning (RingPositionLabel=Para + Locant=Zeta + ring =
Benzene_Tyr → AMBER name "CZ"). Session E's projection functions
take atom indices and return names; they do NOT take ring-position
labels. The chemistry-label-to-atom-name relation is the inverse
projection (used at debug or query time, not at H5 emission).

**Substrate-ready.**

### Gap summary for projection

For ring atoms specifically: zero gaps. The substrate carries the
five-tuple `(Element, Locant, BranchAddress, DiastereotopicIndex,
BackboneRole)` for every ring atom; projection functions render
this to AMBER / IUPAC / BMRB names directly. Pro post-Bundle-C
extension changes ring-position labels but does NOT affect the
identity tuple, so projection stays unaffected.


## Open questions

These are decisions the audit surfaces but does not resolve;
Bundle C drafting (or the user) is the right authority.

### 1. TYM phenolate ring label distinction

The TYM (deprotonated tyrosinate) variant has the same per-ring-atom
labels as TYR (Ipso/Ortho/Meta/Para per benzene); the chemistry
distinction is on the OH atom (PlanarGroupKind switches from
AromaticHydroxyl to AromaticOxide; formal_charge=-1 on OH). The
phenolate's stronger conjugation with the ring π system perturbs
the ring-current intensity (literature variation 5-10 ppm shift
on para-Cζ of phenolate vs phenol; Bovey 1988 ch. 4) — but the
substrate doesn't currently distinguish "ring with neutral OH" vs
"ring with deprotonated O⁻" at the ring-atom level.

**Question.** Should the substrate carry an additional axis on
RingMembership (e.g. `pi_perturbed_by_substituent: bool`, or a
distinct `RingSystemKind::Benzene_Tyr_Phenolate`) to surface the
ring-current chemistry difference? The current TYR vs TYM
distinction is at the variant level (different residue table),
not at the ring level — calculators that read TYR Ring chemistry
do not know whether they're seeing a phenol or phenolate.

**Recommendation deferred.** Bundle C is not the right scope; flag
for a future `Distributed ring current` (planned-audit B4)
calculator extension if quantitative phenolate-vs-phenol
ring-current separation matters.

### 2. TRP perimeter labelling scheme (Gap 4 Option (b))

If Option (b) is chosen for TRP perimeter substrate encoding, what
RingPositionLabel values to use? Two paths:

- **Path 2a — distinct perimeter labels.** `PerimeterIpso`,
  `PerimeterPyrroleBeta`, `PerimeterIndoleNH`, `PerimeterBridgeIn`,
  `PerimeterBridgeOut`, etc. Adds 8-9 new enum values. Per-position
  perimeter chemistry queryable directly.

- **Path 2b — reuse existing labels for secondary.** Each Trp
  atom's secondary `RingMembership.position` reuses a label from
  the existing 14 values (e.g. CG's secondary position = Ipso;
  CD1's secondary = PyrroleBeta; ...). The semantic axis is
  "secondary = perimeter ring", and the position label tells
  perimeter walk role. Saves enum bloat; conflates per-ring labels.

**Recommendation deferred.** Bundle C lands either Option (a)
runtime synthesis or Option (b) substrate extension; if Option
(b), the labelling scheme is the user's call.

### 3. `kAromaticRingTypeCount` named constant placement

Add as `inline constexpr int kAromaticRingTypeCount = 8;` in
`Types.h:185` directly after the enum, OR in a separate `Ring.h`
constant block. Style call.

**Recommendation.** Co-locate with `RingTypeIndex` in `Types.h`
for visibility; calculators including `Types.h` for
`RingTypeIndex::*` see the constant adjacent. The name encodes
"the per-aromatic boundary inside the 9-value RingTypeIndex".

### 4. CreateRing default-fallback discipline

`Ring.cpp:60-61` returns a `PheBenzeneRing` for any unknown
`RingTypeIndex`. Bundle C's substrate-miss FATAL discipline
suggests this should `abort()` instead. Strictly adjacent to
substrate-side scope (factory function lives in `Ring.cpp`, not
in substrate generator); flagged here because Bundle C's plan
explicitly mentions fail-loud-on-substrate-miss in the README §11
of the framing memo. The audit's recommendation:

```cpp
default:
    std::cerr << "FATAL: CreateRing(RingTypeIndex=" << static_cast<int>(type)
              << ") — unknown ring type. Substrate populator emitted an unmapped value.\n";
    std::abort();
```

### 5. Substrate `n_heteroatoms` validation

The substrate generator emits `n_heteroatoms` per `RingMembership`;
spot-check confirmed values match the chemistry (PHE = 0, TYR = 0,
TRP-5 = 1, TRP-6 = 0, HIS = 2, Pro = 1). For HIP, both nitrogens
carry n_heteroatoms = 2 (consistent — both are "ring nitrogens
counted in the imidazole's heteroatom count"). The audit confirms
the field is populated and consistent. No gap.

### 6. Ring atom membership for Pro Cα (chemistry-walk-order N → Cα → ...)

The Pro convention `N → Cα → Cβ → Cγ → Cδ` walks the ring
starting at N, going around to Cδ (which is back to N). The
geometric ring closes via the N-Cδ bond. The convention starts at
N because N is the chemistry-distinct atom (secondary amine,
backbone connector); Cα is next because it's the chiral centre
attached to N AND the next atom in the residue-walk; the rest
walk outward to Cδ which closes the ring.

This convention is OPPOSITE the cyclic walk if read as
"clockwise viewing the ring from one face" — the ring closes via
the implicit N-Cδ bond rather than an explicit cyclic edge. The
SVD ring-normal computation in `Ring::ComputeGeometry` uses the
first three atoms (N, Cα, Cβ) for the cross-product orientation;
the resulting normal direction is well-defined for this walk.

**Status.** Locked decision; not a gap.

### 7. Bundle C scope and per-aromatic-type calculator boundaries

Per the locked Bundle C scope ("substrate-side only; no calculator
updates this pass"), the per-aromatic-type `[8]` arrays in
calculators (BS / HM / PiQuad / Disp + ConformationAtom +
WriteFeatures NPY) stay as-is. Pro Ring's RingTypeIndex=8 falls
outside the `< 8` guard correctly. Pro contributes Intensity=0 to
ring-current calculator output, so calibration column 8 is
identically zero (per the README §9). Any future bundle that
extends per-aromatic-type arrays to per-ring-type will move the
constant from `8` to `kAromaticRingTypeCount` and possibly to
`Count = 9`; the named constant introduced in §1 makes that
change one-edit per site. Substrate-side: no work; flagged so the
naming discipline doesn't bake `8` into the substrate's view of
its calculator clients.

### 8. Pseudoatom membership for Pro ring atoms

Pro ring carbons Cβ/Cγ/Cδ in `kProAtoms` (lines 327-329) carry
`PseudoatomKind::None`; the methylene Hs (HB2/HB3, HG2/HG3,
HD2/HD3) are `PseudoatomKind::Q` with Locant matching. The
substrate is pseudoatom-aware on Pro Hs, which a future
puckering-aware H-position calculator uses (Pro QB/QG/QD). No
gap — consistent with pre-Pro-extension chemistry.


---

End of audit.
