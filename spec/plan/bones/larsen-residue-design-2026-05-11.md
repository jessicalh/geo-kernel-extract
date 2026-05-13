# LarsenResidue — Perception-Driven Tripeptide Model (2026-05-11)

This document supersedes
`spec/plan/typed-tripeptide-topology-design-2026-05-10.md`. The
positional-ordering-table approach in that earlier design is **retired**
before any of its code landed; this document captures what replaces it.

## What changed and why

The 2026-05-10 design proposed:

- Build-time positional ordering tables per `(residue, frame_type)`,
  generated from a hand-authored TOML.
- Runtime composes typed semantics per DFT atom by indexing into the
  positional table.

User concern raised 2026-05-11 (mid-session):

> better check the procs15 source it is possible an old you hallucinated
> things on the way to the db [...] it would probably make sense to
> fully reconstruct that system in some ways. I am suddenly afraid our
> pg data is compromised, unless you can find another way to check it
> [...] which may be implicit in the data [...] what we should really
> have is a LarsenResidue object or something [...] so we are not seeing
> layers of magic order numbers dancing around later.

Two-part concern:

1. **Data trust.** Does the Postgres `tensorcs15` replica preserve the
   atom ordering of the original Larsen ProCS15 Gaussian / ORCA log
   files? If the previous Claude's ingest script perturbed the order,
   our positional-table approach would silently encode the perturbation.
2. **Architectural fragility.** Even if the DB is faithful, encoding
   atom identity as "position 11 is CB, position 12 is BB carbonyl C"
   means every downstream consumer reads the typed slot through an
   opaque positional index. Any future change to the ingest convention
   (new DFT method, different solver, different parser version) breaks
   the magic numbers everywhere.

**Verification done 2026-05-11.** Original Larsen Gaussian log
(`/mnt/expansion/williamsproject/procs15_data/nmrlogs/AAAnmrlog/AAA_4_54_nmr.log`)
inspected at the Standard orientation block. The atom ordering matches
the DB row's `geometry` JSONB byte-for-byte (element sequence
`C O C H H H N H C H C C H H H O …`; atom 11's position is at 1.53 Å
from atom 9's CA, atom 12 is at 1.55 Å from CA and 1.22 Å from atom 16
— confirming atom 11 = CB, atom 12 = BB carbonyl C). The DB is faithful.

The architectural concern stands. **The fix is to derive identity from
chemistry, not from position.**

## Architecture

### LarsenResidue: a first-class peer of AminoAcidType

`src/LarsenResidue.{h,cpp}` is a peer of `src/AminoAcidType.{h,cpp}`,
linked into `nmr_shielding`. It models one residue's-worth of atoms
inside a Larsen tripeptide DFT calculation. Five kinds:

| Kind     | What it represents                                   |
|----------|------------------------------------------------------|
| AceCap   | Acetyl cap: `-C(=O)-CH3` (6 atoms)                   |
| NCapAla  | N-terminal flanking ALA residue (10 atoms)           |
| Central  | Central residue X (any of 20)                        |
| CCapAla  | C-terminal flanking ALA residue (10 atoms)           |
| NmeCap   | N-methylamide cap: `-NH-CH3` (6 atoms)               |

Per-atom record:

```cpp
struct PerAtom {
    AtomMechanicalIdentity identity;     // typed, perception-derived
    int      dft_atom_idx;               // 1-based; for log correlation only
    Element  element;
    Vec3     position;                   // Å, DFT-output frame

    // Pre-decomposed shielding tensor (from DB JSONB; already-irreped at ingest).
    Mat3                  shielding_tensor;
    double                isotropic;
    double                anisotropy;
    std::array<double, 5> t2_components;
};
```

Role-pinned slots (filled by perception, -1 if absent):

```cpp
int N_idx, H_idx, CA_idx, HA_idx, CB_idx, C_idx, O_idx;
```

Predicates / lookups:

```cpp
int  LookupByIdentity(const AtomMechanicalIdentity&) const;
bool HasAllRequiredSlots() const;
```

### LarsenTripeptide: the 5-piece container

```cpp
struct LarsenTripeptide {
    LarsenResidue ace, n_cap, central, c_cap, nme;
    int  TotalAtoms() const;
    std::pair<const LarsenResidue*, int> FindByDftIdx(int) const;
};
```

### Single perception entry point

```cpp
std::optional<LarsenTripeptide> PerceiveLarsenTripeptide(
    const TripeptideDftRecord& rec,
    AminoAcid                  expected_central_residue);
```

Returns `nullopt` if perception fails (e.g., bond graph doesn't
segment into 5 pieces along amide chain; expected central doesn't
match; required slots can't be filled). **Fails loud** — calculator
declines records whose perception fails.

## Perception algorithm

Given `TripeptideDftRecord rec` with `rec.atoms` (positions, elements,
tensors) and `expected_central_residue`:

**Step 1: build bond graph.** Pairwise distance walk over `rec.atoms`.
Element-pair-specific cutoffs from a small lookup table (heavy-heavy
1.85 Å conservative upper bound for C-C single bond + safety; X-H
1.3 Å). The full table:

| Pair         | Cutoff (Å)  | Comment                                |
|--------------|-------------|----------------------------------------|
| C-C, C-N, C-O | 1.85        | Single + double + aromatic            |
| C-S, N-S, O-S | 2.10        | C-S in CYS/MET                         |
| N-H, O-H, S-H | 1.30        | Polar Hs                               |
| C-H          | 1.30        | Aliphatic Hs                           |

Conservatively chosen to capture all canonical bonds without spurious
edges (the typical PM6 / DFT-optimised geometry has bond lengths well
within these cutoffs).

**Step 2: identify amide-backbone chain.** A peptide amide is a
`-C(=O)-N-` motif: an sp2 C with bonds to (one C, one O double-bond,
one N). Walk the graph to find the 4 amides in the tripeptide. Result:
an ordered chain of 5 fragments separated by 4 amide bonds. The
amide bonds are themselves not part of any single piece — they cross
piece boundaries.

**Step 3: segment into 5 pieces.** Each amide-chain endpoint defines
a piece boundary. The 5 pieces are ordered by their position along
the amide chain (the N-end of ACE is the start). Identify ACE by the
methyl-C-bonded-to-carbonyl-C pattern; NME by the methyl-C-bonded-to-
amide-N pattern. Verify size constraints (ACE = 6 atoms, NME = 6,
NCapAla = CCapAla = 10, Central = `AminoAcidType.atoms.size()` for
the expected residue).

**Step 4: per-piece graph isomorphism.** For each piece, run graph
matching against the canonical `AminoAcidType.atoms` inventory (or
the synthetic ACE / NME inventory). Match by element + bond-neighbour
multiset signature. For atoms with degenerate signatures (e.g. ALA's
3 methyl Hs all share the signature "H bonded to one C with 3 H
neighbours"), assign them to the same identity class (DiastereotopicIndex
collapse for methyls; the spatial disambiguation of prochiral Hs is
deferred — see Open Question §1).

**Step 5: emit typed identity per atom.** Run
`ParseAtomName(canonical_atom_name)` (existing function in
`src/generated/LegacyAmberSemanticTables.h`) to produce the typed
`AtomMechanicalIdentity`. The canonical name is the PDB / IUPAC label
from `AminoAcidType.atoms[i].name`.

**Step 6: fill role-pinned slots.** Walk the typed identities and
populate `N_idx`, `H_idx`, `CA_idx`, etc. on each `LarsenResidue`.

**Step 7: verify completeness.** Every canonical atom must have at
least one matched DFT atom. Excess DFT atoms (not in canonical
inventory) → fail loud (the DB row has unexpected content). Missing
canonical atoms → fail loud (the DB row is incomplete for this residue
+ variant; perception declines).

## Integration

### Cache on TripeptideDftRecord

```cpp
struct TripeptideDftRecord {
    // ... existing fields ...

    // Perceived typed model. nullopt if not yet perceived (or perception
    // failed). Populated lazily by QueryNearest on first emit; subsequent
    // calculators read typed fields off the cached LarsenTripeptide.
    std::optional<LarsenTripeptide> larsen;
};
```

`TripeptideDftTable::QueryNearest` calls
`PerceiveLarsenTripeptide(rec, AminoAcidFromLetter(residue_letter))`
before returning the record. Empty / failed perception is logged but
not fatal — the record is still returned with `larsen = nullopt` so
the calculator can decide whether to skip or proceed with a degraded
fallback. (Recommended: skip in production, log + flag in tests.)

### TripeptidePoseAssembler rewrite

`AssembleTripeptide` reads BB N/CA/C indices from the typed slots on
the appropriate `LarsenResidue` piece (`central` for `TripeptidePoseSide::Central`;
`n_cap` for `NTerm`; `c_cap` for `CTerm`). Kabsch input is the typed
positions, not heuristic-identified positions.

Sidechain matching becomes identity-equality:

```cpp
for (auto& dft_atom : larsen.central.atoms) {
    const auto& dft_id = dft_atom.identity;
    const std::size_t protein_ai = LookupByIdentityInProteinResidue(
        protein, residue_idx, dft_id);
    if (protein_ai == kNotFound) {
        ++out.n_substrate_disagreements;
        continue;
    }
    // emit AlignedDftAtom with typed match
}
```

The protein-side typed identity comes from
`LegacyAmberTopology::SemanticAt(protein_ai)`. Both ends now speak
the same vocabulary; matching reduces to a typed-identity equality
check.

The existing `IdentifyCentralBackbone`, `IdentifyAlaCap`,
`IdentifyCTermAlaCap` public APIs are retired (move to anon namespace
or delete). Their failure modes (BB carbonyl C / O confusion in SER,
the O ↔ OG swap) are subsumed by perception.

### CMakeLists.txt

Append `src/LarsenResidue.cpp` to `nmr_shielding` sources. No new
library, no new dependency.

## Test strategy

### Perception correctness (per-residue parity)

`tests/test_larsen_residue_perception.cpp`. For each
`(tripeptide, frame_type)` combination in DB (20 total):

1. Query one row.
2. Run `PerceiveLarsenTripeptide(rec, AminoAcidFromLetter(...))`.
3. Assert `larsen.has_value()` (perception succeeded).
4. Assert all 5 pieces have `HasAllRequiredSlots() == true`.
5. Assert per-piece atom count matches expected.
6. Assert per-canonical-atom-name has at least one match (and
   exactly one for non-prochiral atoms).
7. For BB-bearing pieces (everyone except caps' cap atoms): assert
   N-CA, CA-C, C-O bond distances are in canonical range.

### Source-log forensic test

`tests/test_larsen_residue_against_source_log.cpp`. Parse one
original Gaussian log (`/mnt/expansion/williamsproject/procs15_data/nmrlogs/AAAnmrlog/AAA_4_54_nmr.log`)
directly — read the "Standard orientation" block — and verify that
perception assigns the same typed identity to each atom as it does
for the corresponding DB row. This is the **independent of the DB**
check: if a future DB ingest perturbs the order, this test fails.

### Pose assembler integration

`tests/test_tripeptide_backbone_shielding.cpp` and
`tests/test_tripeptide_neighbor_shielding.cpp` re-run on
`1UBQ_pm6dh3plus.pdb`. Matched-residue count + per-atom residual
distribution must not regress (target: ≥74/76 BB, ≥75/76 neighbor).

### SER OG ↔ O regression

`tests/test_larsen_ser_sidechain.cpp`. Specific check: in a SER
record, after perception, the atom at the canonical OG position
carries `Locant::Gamma` (sidechain hydroxyl O), not
`BackboneRole::CarbonylOxygen`. The 2026-05-10 OG ↔ O swap bug
cannot recur.

## MatchPiece: K=3 Weisfeiler-Lehman + always-relaxed match

Spotted in adversarial review round 3 (codex xhigh, 2026-05-11) as
"Medium: MatchPiece is signature grouping, not graph isomorphism."
Initially marked as known limitation; landed as a fix in the same
session after re-evaluation showed the silent-corruption risk was
real (not just a within-pair swap — the single-round signature class
for PHE/TYR pooled five chemically-distinct atoms CD1/CD2/CE1/CE2/CZ,
allowing cross-locant tensor scramble).

The fix is two-layer:

1. **K=3 Weisfeiler-Lehman signatures in `MatchPiece`** (src/LarsenResidue.cpp).
   Each round's per-atom signature combines the prior round's signature
   with the sorted multiset of neighbours' prior-round signatures
   (FNV-1a 64-bit hashed). After K=3, all chemically-distinct atoms in
   the standard 20 residues are in singleton signature classes. The
   only residual multi-atom classes are graph-automorphic pairs
   (CD1↔CD2, CE1↔CE2, HD1↔HD2, HE1↔HE2 on PHE/TYR; NH1↔NH2,
   HH11/HH12↔HH21/HH22 on ARG; HD21↔HD22, HE21↔HE22 on ASN/GLN;
   methyl Hs that collapse by chemistry anyway).

2. **Per-perceived-atom dispatch in `AssembleCentralTyped`**
   (`src/TripeptidePoseAssembler.cpp`). Round-3 originally landed
   always-relaxed (drop `BranchAddress` + `DiastereotopicIndex` for
   every atom and resolve by nearest-spatial). Round 4 caught that
   always-relaxed was over-broad: it silently dropped CIP-derived
   `BranchAddress` binding for chemistry-distinct branches like ILE
   CG1 vs CG2 (K=3 WL splits them at K=1 — CG1 has a methylene
   extension to CD1, CG2 is a terminal methyl), allowing nearest-
   spatial to swap them under non-canonical chi orientations.

   The current per-atom dispatch (landed 2026-05-11 round 4):

   - `canonical_assignment_ambiguous=false` on `LarsenResidue::PerAtom`
     (singleton WL class — the common case after K=3): STRICT
     identity match. `BranchAddress` and `DiastereotopicIndex` are
     determined by the bond graph per Markley 1998 Fig 1 CIP rules
     and bind to the protein side.
   - `canonical_assignment_ambiguous=true` (multi-atom WL class —
     residual graph-automorphic pairs): RELAXED match dropping
     `BranchAddress` + `DiastereotopicIndex`, with nearest-spatial
     tiebreak within the equivalence class.

   The flag is set in `MatchPiece` when emitting from a canonical
   WL class of size ≥ 2. Mirrored in `scripts/perceive_larsen_tripeptide.py`
   for spec parity (perception output adds
   `canonical_assignment_ambiguous` per atom).

3. **`MatchPiece` returns canonical NODE INDICES, not name strings.**
   The round-3 implementation returned `map<int, std::string>` —
   convenient from the Python POC's shape but a string-as-chemistry-
   carrier crossing the typed-substrate boundary at runtime. Names
   are now constrained to canonical-piece construction only:
   `ParseAtomName` runs once per canonical atom at
   `StampCanonicalIdentities`; `FinalizeAdjacency` then translates
   the name-keyed bond list to index-keyed `adj_by_idx`; runtime WL
   signatures, MatchPiece, and EmitPiece all operate on canonical
   node indices. `EmitPiece` reads typed identity directly from
   `canon.atoms[canon_idx].identity` and dispatches the role-pinned
   slot cache by typed enum (`BackboneRole::Nitrogen` →
   `N_idx`, etc.) — no name comparisons in perception's hot path.
   The user's principle: "Names are only source labels used to
   build canonical topology, not runtime identity and not the
   matching key."

Cost: K=3 WL is ~3000 hash operations per perception call (≈50 atoms ×
~10 neighbour-pairs × 3 rounds × FNV step). Negligible.

Coverage: smoke metrics on 1UBQ_pm6dh3plus.pdb unchanged
(74/76 / 1205/1232 / 0.015 Å mean RMSD), which says all three
dispatch regimes (single-round-strict, round-3 always-relaxed,
round-4 per-atom dispatch) happen to produce the same atom
assignments for 1UBQ's specific geometry. The round-4 fix closes
the silent-swap hole that future MD-frame variance or non-canonical
chi orientations would surface. Tests:

- `test_larsen_residue_perception.cpp` — 20-combo perception parity.
- `test_larsen_residue_against_source_log.cpp` — DB-independent.
- `test_larsen_residue_ser_sidechain.cpp` — SER OG/O distinctness.
- `test_larsen_residue_wl_ambiguity.cpp` (NEW round-4) — ILE
  CG1/CG2 perceive as canonical_assignment_ambiguous=false
  (chemistry-distinct singletons via K=1 WL split); PHE
  CD1/CD2/CE1/CE2 perceive as ambiguous=true (graph-automorphic
  pairs no K splits); PHE CZ perceives as singleton para
  ambiguous=false.

## Round 5 — codex xhigh fourth pass (2026-05-11)

Four findings, all landed.

1. **M1 — Canonical identities now use the generated topology table
   as authority.** Round-4 stamping still ran through `ParseAtomName`
   directly (`src/LarsenResidue.cpp:771` as it was). That can
   diverge from the protein side, where `ComposeAtomSemantic` applies
   a methyl-H pseudoatom collapse (ALA HB1/HB2/HB3 → DI::None) and
   then resolves via `gen::LookupBy(residue, variant_idx, identity)`
   against the substrate generator's output.

   Fix:
   - Split `StampCanonicalIdentities` into `StampCapIdentities`
     (hand-coded for ACE/NME — these are not in the standard-20
     substrate table) and `StampChainIdentitiesViaTable` (for
     `NCapAla`/`CCapAla`/`Central` — the standard-20 chemistry).
   - `StampChainIdentitiesViaTable(p, aa, variant_idx)` mirrors
     `ComposeAtomSemantic` step-for-step: parse name → build
     identity → apply methyl-H collapse via canonical bond graph
     (parent has 3+ H neighbours → clear DI) → `LookupBy(aa,
     variant_idx, ident)` → use the row's identity. Lookup miss is
     FATAL with the same mechanical-identity context the protein
     side emits.
   - Canonical builders thread `(AminoAcid, variant_idx)` through:
     `CanonicalResidue(aa)` uses `kBaseVariantIdx`;
     `CanonicalHisVariant("HID"/"HIE"/"HIP")` resolves to variant
     index 0/1/2 (matching `Residue::protonation_variant_index`
     convention and `AminoAcidType` variant ordering).
   - `FinalizeAdjacency(p)` is now called BEFORE the stamper for
     chain pieces so `adj_by_idx` is available for the methyl-H
     parent-degree count.
   - Removed unused `CanonicalIdentity()` helper.

2. **M2 — Python POC framing.** The header comments in
   `src/LarsenResidue.h:21` and `src/LarsenResidue.cpp:1` previously
   called `scripts/perceive_larsen_tripeptide.py` "the validated
   spec for this C++ port". The Python script still returns
   canonical name strings from `match_piece` (predates the typed-
   substrate work) and is no longer the spec — the C++ object model
   is the runtime authority. Comments updated to call the Python
   script the original prototype that validated the bond-graph + WL
   algorithm, useful as a diagnostic tool, not a normative spec.

3. **L1 — Neighbor test NaN handling.**
   `tests/test_tripeptide_neighbor_shielding.cpp` previously read
   `.norm()` unconditionally and compared against 3.0; under the
   M4 NaN-fill contract, `Vec3{NaN}.norm()` is NaN and
   `NaN > 3.0` is false — so a direction being absent everywhere
   would silently pass the extreme-bound check. Fixed: count
   `n_finite_prev` / `n_finite_next` separately, only fold finite
   vectors into the max/extreme tallies, and assert non-trivial
   finite coverage in each direction (> 100 finite atoms).

4. **L2 — Assembler header docstring tightened.**
   `src/TripeptidePoseAssembler.h` previously described
   `validation_threshold_A` as "atoms whose residual_distance
   exceeds this are excluded from aligned_atoms". After round 1,
   the central path only records the threshold diagnostically
   (does NOT reject — residual_vec is the load-bearing ML feature
   per `feedback_residual_as_ml_feature`); only the cap path still
   gates. The docstring now separates the two regimes explicitly so
   future work doesn't resurrect the old rejection model on the
   central path.

Smoke unchanged on 1UBQ_pm6dh3plus.pdb after round 5: BB 74/76 /
1205/1232 atoms (97.8%) / mean RMSD 0.015 Å / max 0.048 Å;
Neighbor 76/76 / 522 accumulations. 217 non-tripeptide + 7
tripeptide structure tests pass; 95 Python SDK tests pass.

## Round 6 — deferred-items resolution sweep (2026-05-11)

Five items surfaced after the codex R5 review showed clean — all
landed in this round.

1. **DB single-variant verification + HIS warn.** `psql` query
   against `raw_dft_calculations` confirmed each tripeptide letter
   has uniform `n_atoms`: non-HIS titratable residues (CYS, ASP,
   GLU, LYS) all carry the default AMBER form, and AHA carries 50
   atoms = 32 caps + 18 central = HIP only. The
   `kBaseVariantIdx` assumption in `StampChainIdentitiesViaTable`
   is verified for non-HIS; for HIS, `TripeptideBackboneShieldingResult::Compute`
   and `TripeptideNeighborShieldingResult::Compute` now warn loudly
   when a protein-side HIS has `protonation_variant_index ∈ {0, 1}`
   (HID/HIE) because perception will fail and σ_BB^i / Δσ_BB will
   be absent for that residue. The trigger to revisit the warn is
   future tensorcs15 ingest gaining HID/HIE rows.

2. **FATAL message anchor.** `StampChainIdentitiesViaTable`'s
   FATAL `fprintf` now names `LarsenResiduePerceptionTest.AllCombinationsPerceiveCleanly`
   as the positive-coverage anchor — every chain atom in every
   standard residue × variant is exercised through `LookupBy`,
   so a passing test suite means the FATAL is unreachable from the
   standard-20 substrate. The path can only fire on
   development-time edits that diverge one side without the other.

3. **`ChiFallbackIsDeterministic` test.**
   `tests/test_larsen_residue_wl_ambiguity.cpp` gains a regression
   test verifying that two consecutive `QueryNearest('R', -180,
   -180, ..., n_chi_axes=0)` calls return the same calc_id —
   confirming round-4 M3's `ORDER BY calc_id ASC` makes
   chi-fallback row selection deterministic.

4. **HBondHα memory: sidechain-O acceptor trigger.**
   `project_hbond_halpha_design` memory entry rewritten — the
   "deferred unless extended later" language replaced with
   explicit Phase 1 scoping + post-calibration trigger condition.

5. **Phase-2 prochiral methylene H disambiguation: trigger
   condition.** Open Question §1 in this doc rewritten — Phase-2
   stays unimplemented; the explicit revisit trigger is
   post-Stage-2 calibration showing systematic pro-R/pro-S bias on
   methylene Hα/HB/HG/HD/HE pairs (which would indicate spatial
   tiebreak is mis-assigning the locant under MD chi variance).

Smoke unchanged: BB 74/76 / 1205/1232 / 0.015 Å mean; Neighbor
76/76 / 522. 217 non-tripeptide + 8 tripeptide structure tests
pass; 95 Python SDK tests pass.

## Open questions

### 1. Prochiral methylene H disambiguation

ALA's HB1/HB2/HB3 are a methyl (DiastereotopicIndex collapses; any
of the three can match either of the three protein-side Hs).
For prochiral methylenes (HB2/HB3 on ARG/SER/GLU/etc., HD2/HD3,
HE2/HE3, HZ2/HZ3, etc.), the IUPAC convention assigns position-2
(pro-R) and position-3 (pro-S) by CIP priority, which requires
spatial reasoning.

**Phase-1 decision (LANDED via round-4 ambiguity flag + R5
substrate-table grounding):** perception puts prochiral methylene
Hs (HB2/HB3 on ARG/SER/etc., HD2/HD3, HE2/HE3, HZ2/HZ3) in a
multi-atom K=3 WL signature class. `canonical_assignment_ambiguous=true`
on each. `StampChainIdentitiesViaTable` then looks up the canonical
row in the generated substrate table — the per-DFT-atom typed
identity reflects what the table says about each canonical name
(`DiastereotopicIndex::Position2/Position3` AND
`ProchiralStereo::ProS/ProR` are both encoded in the generated SER
row at `src/generated/LegacyAmberSemanticTables.cpp:349-350`,
verified 2026-05-11). At match time, `AssembleCentralTyped` uses
relaxed identity (drops `BranchAddress` + `DiastereotopicIndex`)
plus nearest-spatial tiebreak — the Kabsch-rotated DFT positions
land near the protein's HB2/HB3 in the same pro-R/pro-S
arrangement when chi geometry is faithful, so spatial tiebreak
picks the correct locant.

**Phase-2 status: NOT IMPLEMENTED, NOT BLOCKED. Trigger:**
post-Stage-2 calibration, when per-atom-type stratification of the
residual against DFT shows a systematic pro-R vs pro-S asymmetry on
methylene Hα/HB/HG/HD/HE pairs. The asymmetry would indicate that
spatial tiebreak is mis-assigning the locant at some non-trivial
fraction of frames — most likely under MD frame variance where the
chi angle differs enough from Larsen's 20° grid resolution to flip
the physical pro-R/pro-S H positions relative to the DFT
canonical assignment. Until that asymmetry is measured, the
spatial tiebreak is the right behaviour and Phase-2 stays
unimplemented.

**Phase-2 work (when triggered):** add bond-graph CIP perception in
`LarsenResidue.cpp` that assigns `ProchiralStereo::ProR/ProS` per
methylene H using neighbour priority walking, mirroring RDKit's
CIPLabeler. Stamp the result onto `LarsenResidue::PerAtom`. The
generated table already carries the canonical pro-R/pro-S labels;
the perception side just needs to discriminate two graph-isomorphic
Hs by spatial priority of the heavy-atom branches around their
parent carbon.

### 2. PRO ring perception

PRO closes a ring via N-CD bond. Perception must recognise the ring
and not double-count atoms. The N atom is shared between backbone
and the pyrrolidine ring; per the existing `RingSystemKind::Pyrrolidine_Pro`
substrate field this is encoded.

**Decision:** PRO perception uses the standard amide-segmentation
+ graph-iso path; the ring closure is detected as a cycle in the
bond graph and tagged via the existing `LegacyAmberTopology` PRO
substrate. No special handling needed in `LarsenResidue` itself —
the canonical PRO atoms are matched normally.

### 3. ACE / NME identity assignment

ACE atoms (carbonyl-C, carbonyl-O, methyl-C, 3×methyl-H) and NME
atoms (N, amide-H, methyl-C, 3×methyl-H) are not in the existing
`LegacyAmberSemanticTables.cpp` standard-20 residue list. Their
typed identity needs synthesized values:

- ACE carbonyl-C: `BackboneRole::CarbonylCarbon`, `Locant::None`.
- ACE carbonyl-O: `BackboneRole::CarbonylOxygen`, `Locant::None`.
- ACE methyl-C: `BackboneRole::None`, `Locant::None`, with
  `PseudoatomKind::M` on the methyl-Hs.
- NME amide-N: `BackboneRole::Nitrogen`, `Locant::None`.
- NME amide-H: `BackboneRole::AmideHydrogen`, `Locant::None`.
- NME methyl-C and methyl-Hs: as for ACE.

These are hand-coded in `LarsenResidue.cpp` (the 5 kinds need only
~12 explicit identity assignments total for ACE+NME). No round-trip
through the generator.

### 4. Failure logging

When perception fails for a record, what gets logged + what does
the calculator do?

**Decision:** structured warn via `OperationLog` with the record's
`calc_id`, `tripeptide`, `frame_type`, and the perception-failure
reason. The calculator skips the residue and continues. The user
can audit per-residue skip rates; if a particular residue + frame
type fails consistently, that's a perception bug to fix, not a
data bug.

## File layout

**New:**

| Path                                            | Lines  | Purpose                                    |
|-------------------------------------------------|--------|--------------------------------------------|
| `src/LarsenResidue.h`                           | ~120   | LarsenResidue, LarsenTripeptide, Perceive entry |
| `src/LarsenResidue.cpp`                         | ~400   | Perception algorithm; hand-coded ACE/NME identities |
| `tests/test_larsen_residue_perception.cpp`      | ~120   | Per-(residue, frame_type) DB query + perception parity |
| `tests/test_larsen_residue_against_source_log.cpp` | ~80 | Independent of DB; one original Gaussian log |
| `spec/plan/larsen-residue-design-2026-05-11.md` | this   | Durable design (this doc)                  |

**Modified:**

| Path                                  | Change                                                |
|---------------------------------------|-------------------------------------------------------|
| `src/TripeptideDftTable.h`            | Add `optional<LarsenTripeptide> larsen;` on Record    |
| `src/TripeptideDftTable.cpp`          | Call `PerceiveLarsenTripeptide` in `QueryNearest`     |
| `src/TripeptidePoseAssembler.{h,cpp}` | Rewrite to read typed slots from LarsenResidue        |
| `tests/test_tripeptide_backbone_shielding.cpp` | Updated expectations (no behavior change expected)  |
| `tests/test_tripeptide_neighbor_shielding.cpp` | Updated expectations (no behavior change expected)  |
| `CMakeLists.txt`                      | `+ src/LarsenResidue.cpp`                             |

**Retired (no code lands):**

- `data/topology/tripeptide_orderings.toml` — never authored
- `data/topology/governing_dihedrals.toml` — never authored
- `src/generated/TripeptideOrderingTables.{h,cpp}` — never generated
- The Phase-1 sub-tasks of the old design — deleted from task list

The dumper script (`scripts/dump_tripeptide_orderings.py`) is **kept**
as a forensic diagnostic: it can compare DB ordering against a fresh
re-parse of source logs to detect drift over time. Not used as input
to generator code.

## Why this is the right architectural shape

1. **Identity is chemistry, not position.** Bond graph + canonical
   `AminoAcidType` are the chemistry; positions are observation. We
   match observation to chemistry, not chemistry to position-in-array.

2. **DB ingest can change without breaking downstream.** If a future
   DB ingest reorders atoms or adds new fields, perception still
   produces the same typed identity because the chemistry doesn't
   change. The positional ordering becomes a side-channel diagnostic.

3. **Substrate-side discipline unchanged.** The protein side already
   uses typed identity from `LegacyAmberTopology::SemanticAt(ai)`;
   this design extends the discipline to the DFT side. Cross-substrate
   matching is now typed-identity ↔ typed-identity (per
   `feedback_two_path_validation`).

4. **No build-time generator dependency on chemistry strings.** The
   existing `LegacyAmberSemanticTables` generator at
   `tools/topology/build_semantic_tables` continues to be the SINGLE
   place chemistry strings cross the build-time boundary. Runtime
   `LarsenResidue` consumes typed atoms from that table; perception
   is a runtime layer that produces typed atoms from raw DFT data.

## Implementation order (this session, if time permits)

1. ✅ Write this design doc.
2. Python POC: `scripts/perceive_larsen_tripeptide.py`. Run on all
   20 DB combos + sanity-check against original log for AAA. Validates
   the algorithm before C++ port.
3. `src/LarsenResidue.h`.
4. `src/LarsenResidue.cpp`.
5. Wire into `TripeptideDftRecord` + `TripeptideDftTable::QueryNearest`.
6. Rewrite `TripeptidePoseAssembler::AssembleTripeptide`.
7. Tests: perception parity, source-log forensic, pose assembler
   integration, SER regression, BB + Neighbor smoke.
8. Update CLAUDE.md + memory.

Sizing: 4–6 hours, possibly two sessions. Phase 1 (design doc + Python
POC + C++ scaffolding) likely fits this session; Phase 2 (pose
assembler rewrite + tests + smoke) likely next session.

## Provenance

- 2026-05-10 design (`typed-tripeptide-topology-design-2026-05-10.md`)
  produced by Plan agent + user pose-metadata enrichment. Retired
  before any of its code landed.
- This document produced 2026-05-11 after user raised data-trust +
  architectural-fragility concerns mid-session. DB ordering verified
  against original Larsen Gaussian log; pivot to perception captured
  here.
- The Python POC + C++ implementation that follow consume this
  document as authoritative.
