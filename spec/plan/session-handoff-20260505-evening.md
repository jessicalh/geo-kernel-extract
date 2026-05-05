# Session handoff — 2026-05-05 evening

**Pickup brief.** This evening session built on the Session B morning
handoff (`spec/plan/session-handoff-20260505.md`). It produced the
consolidated chemistry research the user explicitly asked for ("one
deep research pass across the full set, then proceed carefully"),
took the resulting reference document through critical review by an
empowered LLM agent, captured the application-time decision points
in a dependencies file, and applied the agreed corrections + enum
extensions in one coherent commit.

The encoding of the 19 remaining residues + 10 variants is the next
session's job. It is now mechanical transcription against signed-off
chemistry inputs, not chemistry research.

The user requested a load-bearing handoff (their flagged risk:
direction-and-perception drift across context-to-context transitions).
This document is that.

---

## State at end of session

**HEAD:** `ab3bd5e` on `master`. The morning handoff commit `b079d4a`
plus three evening commits (`bb4584d`, `56fb6b9`, `ab3bd5e`) are
local-only pending push (matches the prior session's state — push is
on the user's schedule, not autonomous).

**This session's three commits, oldest first:**

```text
bb4584d  Phase 5 substrate: residue reference + critical review +
         dependencies + enum extensions
56fb6b9  Phase 5 substrate: encode standard 20 residues; variants pending
ab3bd5e  Phase 5 substrate: emit AMBER variant tables via parent-CCD +
         delta architecture
```

**Tests:** ctest 352/352 pass (347 prior baseline + 5 StringBarrier).
Generator regenerates byte-identical PHE output. Working tree clean.

---

## What this session produced

### Three plan documents (the encoder's source materials)

1. **`spec/plan/topology-residue-reference-2026-05-05.md`** (898 lines)
   Per-(residue, variant, atom) chemistry encoding for all 20 standard
   residues + 10 AMBER protonation variants + 4 terminal states. PHE
   rows match `SynthesisedForPhe` exactly. The "what" of encoding.

2. **`spec/plan/topology-reference-review-2026-05-05.md`** (207 lines)
   Critical review by an empowered LLM agent: 7 ACCEPTs, 1 ACCEPT-WITH-FLAG,
   1 REJECT (Asn/Gln amide H E/Z inversion — caught real bug), 4 MODIFY
   (LYN HZ atom names, AmineNH addition, AromaticOxide addition, terminal-H
   pseudoatom doc). The "agreed signoff" that the encoder operates against.

3. **`spec/plan/topology-encoding-dependencies-2026-05-05.md`** (new)
   Application-time decision manifest. Six sections:
   - §A vocabulary preconditions (already applied this session)
   - §B external-authority dependencies (encoder MUST verify)
   - §C reading conventions consumers MUST know
   - §D locked architectural decisions (do not relitigate)
   - §E reference-doc corrections (already applied this session)
   - §F ARN remains a flag, not a verdict

### Two enum values added to `src/SemanticEnums.h`

- `PolarHKind::AmineNH = 11` — neutral primary-amine NH (LYN HZ2/HZ3,
  NTERM_NEUTRAL H1/H2). Distinguished from `AmmoniumNH` per pKa,
  exchange behaviour, and SHIFTX2 separate-model evidence.
- `PlanarGroupKind::AromaticOxide = 9` — TYM Oη after deprotonation.
  Phenolate -O- conjugates more strongly with the ring π system than
  neutral phenol-OH does.

`OtherPolarH` renumbered from 11 → 12 (transparent to consumers since
runtime tables are emitted via the literal-emitter switch cases, all
of which were extended in the same commit).

### Three documentation comment blocks

In `src/SemanticEnums.h`:
- `RingPosition` struct: orthogonality of `Locant` and ring membership
  (Pro Cα = the canonical example).
- `RingPositionLabel::Ortho1`/`Ortho2`/`Meta1`/`Meta2`: extension for
  Trp 6-ring perimeter convention (synthesised; Markley does not
  publish ipso/ortho/meta/para for indole 6-ring).
- `PseudoatomKind`: terminus convention (`Q` with `Locant::None` for
  N-terminal H1/H2/H3) — staying within IUPAC M/Q/R taxonomy as a
  discipline call.

### One dead-code removal in `tools/topology/build_semantic_tables.cpp`

`SynthesisedForPhe` line 794 had `|| atom_id == "HG"` matching a
non-existent PHE atom — flagged by reviewer as forward-compat scaffolding
that risks future TYR generalisation accidents. Removed; PHE output
unchanged (verified byte-identical regeneration).

### Reference-document corrections (11 mechanical edits)

Per dependencies §E:
- Asn HD21 Z→E, HD22 E→Z (CIP priority on Cγ is OD1 > CB).
- Gln HE21 Z→E, HE22 E→Z (same logic on Cδ).
- LYN: invert "HZ3 removed" → "HZ1 removed" (per ff14SB).
- LYN HZ2/HZ3, NTERM_NEUTRAL H1/H2: AmmoniumNH → AmineNH.
- TYM Oη: AromaticHydroxyl → AromaticOxide.
- §2 charge-distribution conventions: rewrote Asn/Gln Z/E sub-clause to
  explain residue-specific CIP-priority origin (Asn/Gln "1=E"; Arg
  "1=Z" because of priority on the planar carbon).
- Disagreement log entries 1, 2, 8, 9 marked RESOLVED with action.

---

## Update: encoding completed in this session via supervised agent

After the original handoff was drafted, the user proposed using an
encoding agent under supervisor oversight rather than saving state for
a fresh session. Two agent runs completed in this session:

**Run 1 — encoded standard 20 + per-variant patch functions** but iterated
variant codes against the CCD assuming the codes had matching CCD
entries. They do not — the codes (HID/HIE/HIP/ASH/GLH/CYX/CYM/LYN/ARN/
TYM) coincide with unrelated small molecules in CCD by historical
accident. The 10 variant tables emitted contained garbage data
encoded against the wrong molecules.

**Strip + standard-20 commit (`56fb6b9`)** removed the variant entries
from the iteration loop and dispatch; per-variant `SynthesisedFor<Variant>`
functions kept as `[[maybe_unused]]` ready for the corrected
architecture; `topology-encoding-dependencies-2026-05-05.md` §G added
to capture the architectural lesson.

**Run 2 — variant tables emitted via parent-CCD + delta architecture**
(commit `ab3bd5e`). The agent discovered the CCD parent entries
already contain protonated/imidazolium variants, so most variants are
atom-removal cases (HIE removes HD1; HID removes HE2; CYX/CYM remove
HG; LYN removes HZ1; ARN removes HE; TYM removes HH) plus field-deltas
via the existing `SynthesisedFor<Variant>` patches with a new
`charge_override` field on `SynthesisedFields`. HIP/ASH/GLH have no
atom-removals (parent CCD atoms = variant atoms; just field changes).

Generated table: 30 residue tables (20 standard + 10 variants).
Standard 20 byte-identical. PHE byte-identical. ctest 352/352.

**Phase 5 substrate is COMPLETE.** Next session's work is integration
+ test coverage, not encoding.

---

## What's next (sequenced)

The encoding work is done. Next is integration + test coverage.

### Next-session priority 1: Coverage test (Task #10)

`tests/test_legacy_amber_semantic_tables.cpp` (or similar). Asserts
programmatic invariants on the 30 emitted tables. Suggested asserts:

- **Standard 20**: each has expected backbone atoms (N, CA, C, O,
  plus H for non-Pro); each non-Pro backbone H has
  `PolarHKind::BackboneAmide`.
- **Variant deltas** exercise their critical fields:
  - HID HD1 = `ImidazoleNH`; Nδ1 = `Heteroatom_NH`; Nε2 = `Heteroatom_NoH`.
  - HIE HE2 = `ImidazoleNH`; Nε2 = `Heteroatom_NH`; Nδ1 = `Heteroatom_NoH`.
  - HIP both Hδ1 + Hε2 = `ImidazoleNH`; Nε2 formal_chg=+1.
  - LYN: 24 atoms; HZ1 absent; HZ2/HZ3 = `AmineNH`; Nζ formal_chg=0.
  - LYS: HZ1/HZ2/HZ3 = `AmmoniumNH`; Nζ formal_chg=+1.
  - TYM: Oη = `AromaticOxide`, formal_chg=-1; HH absent.
  - ASH HD2 = `CarboxylOH`; Oδ2 formal_chg=0.
  - GLH HE2 = `CarboxylOH`; Oε2 formal_chg=0.
  - CYM: HG absent; Sγ formal_chg=-1.
  - CYX: HG absent.
  - GLY: HA2 = `ProR`; HA3 = `ProS` (the inversion).
  - TRP: Cδ2 + Cε2 both have primary AND secondary RingMembership.
- **Markley alternation**: β-onward methylenes have HB3=ProR,
  HB2=ProS (sample LYS, ARG, MET, ASP).

Tests can dispatch via a small lookup helper or directly index the
generated `kFooAtoms` / `kFooAtoms_VARIANT` arrays.

### Next-session priority 2: Runtime integration

Add `src/generated/LegacyAmberSemanticTables.cpp` to the
`nmr_shielding` source list in the main `CMakeLists.txt`, plus a
small lookup function `LookupAtomSemantic(AminoAcid, atom_local_idx,
variant_idx) -> const AtomSemanticTable*`.

### Then: Session C+D paired (substrate extension + ConformationResult)

Per the morning handoff. Extend `LegacyAmberTopology` with the 14
typed fields populated via the lookup at `Protein::FinalizeConstruction`.
New `PlanarGeometryResult` ConformationResult for per-frame
ω, ring-flip, sp2 pyramidalisation, ring-pucker geometry.

### Then: Session E (IUPAC + BMRB projections)

User identified as "a must" for analysis to continue. Pure functions
on `LegacyAmberTopology`.

---

## File map

```text
spec/plan/topology-fields-research-2026-05-05.md          [original research dossier (Session B morning)]
spec/plan/topology-substrate-implementation-plan-2026-05-05.md  [the central spec]
spec/plan/topology-residue-reference-2026-05-05.md        [THIS session — encoder source]
spec/plan/topology-reference-review-2026-05-05.md         [THIS session — empowered review]
spec/plan/topology-encoding-dependencies-2026-05-05.md    [THIS session — application-time deps]
spec/plan/session-handoff-20260505.md                      [Session B morning handoff]
spec/plan/session-handoff-20260505-evening.md              [this file]

src/SemanticEnums.h                                        [extended this session: AmineNH, AromaticOxide, doc comments]
src/generated/LegacyAmberSemanticTables.cpp                [PHE only; unchanged this session]
src/generated/LegacyAmberSemanticTables.log.txt            [PHE only; unchanged this session]

tools/topology/build_semantic_tables.cpp                   [dead HG-clause removed; literal-emitter switches extended]

tests/test_string_barrier.cpp                              [unchanged this session]
```

---

## Pitfalls / concerns for the encoder (next session)

1. **The reference doc is markdown, not machine-readable.** Encoding
   means hand-translating tables to C++ switch cases. Error-prone;
   coverage test (Priority 2) is the safety net. Land coverage early,
   not after all 19 residues — partial-coverage failures should fire
   per-residue-batch.

2. **`SynthesisedFor<Residue>` only populates the synthesised +
   Markley-derived fields**, not the mechanical (Locant, BranchAddress,
   DiastereotopicIndex) or algorithmic (aromatic, equivalence_class)
   fields. The latter come from the atom-name parser and RDKit at
   build time. `BuildAtomSemanticEntry` (lines 888+) is where the
   three sources combine.

3. **HIE-default-for-HIS is invisible in code** but locked in
   dependencies §D.1. The HIS rows in the reference ARE the HIE rows;
   the encoding must respect this. If a future session changes the
   AMBER default, the entire HIS family encoding regenerates.

4. **ARN remains a flag, not a verdict** (dependencies §B.1, §F).
   Encoder commits "HE removed" UNLESS `AmberAminoAcidVariantTable["ARN"]`
   in `src/AminoAcidType.cpp` says otherwise. Check at code-gen time.

5. **Charge placement is Lewis-localised** (dependencies §C.1). Asp
   Oδ2 = -1, Oδ1 = 0 (not -0.5 each). HIP Nε2 = +1, others = 0. Arg
   Nε = +1, others = 0. Calculators wanting symmetrised charges must
   transform downstream.

6. **The Markley alternation rule for prochiral Hs is HAND-ENCODED in
   the reference**, not derived from RDKit. RDKit's CIPLabeler labels
   chiral centres (CA → S for L-amino acids), NOT prochiral Hs. The
   encoder reads the reference's `prochiral` column directly per atom;
   does NOT attempt to derive it algorithmically.

7. **Trp 6-ring perimeter labels are synthesised, not literature-canonical**
   (dependencies §C.3). If a downstream calculator stratifies Trp ring
   atoms by ring position with a different convention, alignment is
   needed. No literature contradiction; just non-canonical.

8. **The reviewer's verdict is signed off but enforcement is on the
   encoder.** No automated check ties C++ encodings back to reference
   markdown. The discipline is: encoder reads dependencies + reference
   each time, cross-checks against the reference's quality-bar
   self-audit (Section 7 of the reference).

9. **clangd LSP false positives** persist on `build_semantic_tables.cpp`
   for RDKit headers (editor scans the runtime build, not the opt-in
   generator build). Ignore. Documented in
   `feedback_test_invocation_via_ctest` memory.

10. **Don't investigate the reverted IUPAC-topology episode** (memory
    `feedback_no_autonomous_file_ops` + CLAUDE.md preamble). Substrate
    work is the architectural successor; the prior attempt's code is
    not relevant.

---

## How to verify state in a fresh session

```bash
# 1. HEAD matches.
git log --oneline -5
# Expect: ab3bd5e, 56fb6b9, bb4584d, b079d4a, 721e681

# 2. ctest baseline.
ctest --test-dir build -j 8 | grep -E "tests passed|tests failed"
# Expect: 100% tests passed, 0 tests failed out of 352

# 3. String barrier intact.
ctest --test-dir build -R StringBarrier
# Expect: 5/5 pass.

# 4. Linker barrier still holds.
nm build/libnmr_shielding.a | grep -c '_ZN5RDKit'
# Expect: 0.

# 5. Generated tables emit 30 residues.
grep -c "^// === " src/generated/LegacyAmberSemanticTables.cpp
# Expect: 30 (20 standard + 10 variants).

# 6. Generator regenerates byte-identical for the standard 20.
./build-gen/tools/topology/build_semantic_tables \
    --ccd data/ccd/components.cif \
    --output /tmp/regen.cpp \
    --log /tmp/regen.log
diff /tmp/regen.cpp src/generated/LegacyAmberSemanticTables.cpp
# Expect: no output (byte-identical).

# 7. New enum values present.
grep -E "AmineNH|AromaticOxide" src/SemanticEnums.h
# Expect: 4+ hits (1 enum line + 1 docstring per value).

# 8. Five plan docs in place.
ls spec/plan/topology-*-2026-05-05.md spec/plan/session-handoff-20260505*.md
# Expect: 5 topology docs + 2 handoffs (morning + evening).
```

---

## Closing note

The user's pacing concern is real. Three sessions have now landed on
this work (Session B morning: architecture + PHE proof; this session:
research + review + signoff + corrections; next session: encoding).
The architecture is settled, the chemistry is signed off, the
dependencies are captured. The encoding is data-entry against those
inputs, gated by the coverage test.

No named debt carries forward. The session can pick up cold by reading,
in this order:

1. This file (handoff).
2. `spec/plan/topology-encoding-dependencies-2026-05-05.md` (what to honour).
3. `spec/plan/topology-residue-reference-2026-05-05.md` (the data).
4. `spec/plan/topology-reference-review-2026-05-05.md` (only if a verdict is questioned).

If Section 7 of the reference (quality-bar self-audit) doesn't match
what the encoder produces, that's a flag — fix the encoding, not the
reference. The reference is the agreed truth.
