# Session handoff — 2026-05-05

**Pickup brief.** This session opened with the user asking for a
docs-and-direction triage of Phase 5 (the LegacyAmberTopology
semantic-substrate work) and ended with the architecture validated
end-to-end on one residue (PHE), the string barrier enforced at the
linker level, and the generator emitting a typed-enum table that
compiles standalone. Five commits landed; ctest + pytest stay at
their pre-session baselines (347/347 + 64/68).

The user asked explicitly that this handoff be **load-bearing for
context-to-context continuity**, since they are managing the risk of
direction-and-perception drift across sessions. This document is
structured to support that: explicit "do not relitigate" markers,
the architectural rationale captured alongside the decisions, and
the file map so a fresh session knows where things live.

---

## State at end of session

**HEAD:** `721e681` on `master`. `origin/master` matches the local
master *up to the prior session* (`b9636f0`); these five commits are
local-only pending push. Working tree clean.

**This session's commits, oldest first:**

```text
7390a6d  Quick wins: EDR path + aimnet2_charge_sensitivity retire
d915fb4  Phase 5 scaffolding: SemanticEnums + topology generator
6d4aa12  Topology generator: cifpp + RDKit pipeline (PHE round-trip)
b3e281b  Topology generator: per-field population for PHE
721e681  Topology generator: emitter + string-barrier test
```

**Tests at end of session.** All five `StringBarrier.*` tests pass.
JobSpec suite re-verified (24/24). The full ctest run was not
re-executed at session end because the substrate work has not yet
been integrated into the runtime build path (the generated table is
committed but no runtime code consumes it yet) — so existing tests
are unaffected by the generator.

---

## Architecture as decided (do not relitigate)

The user pushed back hard during this session against my reflexive
splits-and-second-guesses. The following decisions are LOCKED for
continuity. A future session can revisit them but only after
explicit user decision; do not unilaterally re-open.

### A. The substrate carries 14 typed semantic fields

Per `OBJECT_MODEL.md:2486-2498` plus the bonus fields surfaced in
`spec/plan/topology-fields-research-2026-05-05.md`:

```text
Locant                  -- Greek-letter sidechain position
BranchAddress           -- two-level (outer, inner) for Arg eta-Hs
DiastereotopicIndex     -- IUPAC numeric label 2/3 on prochiral methylenes
ProchiralStereo         -- CIP R/S
PlanarGroupKind         -- 8-value classification of sp2 / planar groups
PlanarStereo            -- E/Z canonical label
PseudoatomMembership    -- Markley M/Q/R + locant + branch + super-group flag
PolarHKind              -- 12-value classification of exchangeable Hs
RingPosition            -- structured (primary, secondary) of RingMembership
aromatic                -- bool
formal_charge           -- int8_t
is_exchangeable         -- bool, derived from PolarHKind
equivalence_class       -- uint8_t, RDKit canonical rank
```

Hybridisation and BondOrderMask are computed by the generator and
recorded in the log but are NOT in the runtime AtomSemanticTable
record. Add them to the runtime when a calculator demands them.

### B. The string barrier is enforced at the linker level

`libnmr_shielding.a` does not link RDKit or gemmi. The only
chemistry-string library reachable from runtime is cifpp at the
existing PDB-loading boundary (`PdbFileReader.cpp`,
`DsspResult.cpp`). The discipline is structural: contributors
cannot accidentally call RDKit from runtime because the symbols are
not present in the archive.

Layered enforcement (each tightens the barrier without depending on
the others):

```text
Layer 1 (linker)    -- nm of libnmr_shielding.a shows zero RDKit symbols
Layer 2 (header)    -- zero src/ headers #include cifpp/RDKit/gemmi
Layer 3 (module)    -- cifpp #include + cif:: namespace use confined
                       to PdbFileReader.cpp + DsspResult.cpp; RDKit
                       and gemmi banned everywhere in src/
Layer 4 (generated) -- src/generated/*.cpp contains no std::string
                       literals (only typed-enum identifier text)
```

`tests/test_string_barrier.cpp` enforces all four layers.

### C. RDKit is the citable primary for algorithmic fields

Methodology decision (the user's "RDKit is citable even if we are
wrong, honestly"): a citable, version-pinned, reproducible
algorithmic source wins as primary over hand-curated tables and
historical text. The thesis methods chapter cites
"RDKit version 2023.09.6 CIPLabeler" as a defensible source.
Markley 1998 and BMRB are **substantive** cross-checks --
disagreements get typed-logged in the witness stash for downstream
stats analysis to query.

Per-field source precedence is captured in
`spec/plan/topology-substrate-implementation-plan-2026-05-05.md`
"Reconciliation precedence" table.

### D. Synthesised fields are first-class chemistry contributions

The user explicitly corrected my "weaker provenance" framing. The
synthesised fields (`PlanarGroupKind`, `PolarHKind`,
`PseudoatomMembership`, `RingPosition` labels for 5-ring/fused-ring
positions) are:

- Grounded in cited chemistry literature (Pauling 1951, Cantor &
  Schimmel 1980, Joule & Mills 2010, Vollhardt & Schore 2018,
  Wuethrich 1986, Englander 2008).
- Encoded with per-enum-value docstring citations in
  `src/SemanticEnums.h`, same standard as `Ring::Intensity()`'s
  Case 1995 citation in `src/Ring.cpp`.
- Owned by us. The taxonomy is defensible because the chemistry is
  cited.

Do NOT pre-emptively diminish this work. The ring-currents work in
`src/Ring.cpp` is the precedent: hand-built typed properties with
literature citations are first-class scientific contributions, not
fallbacks.

### E. The substrate-vs-conformation split for planar stereo

The user said "we need conformation, conformation is not just binary
but angle." The literature (SHIFTX2: ω accounts for 5-14% of
backbone shielding signal, treated as continuous; ProCS15: locks
ω=180 in its substrate scan and is wrong to do so) supports:

- **Substrate side** (this work): `PlanarGroupKind` enum (which atoms
  are in a planar group) + `PlanarStereo` E/Z canonical label.
- **Conformation side** (next session, NOT in this work): per-frame
  `PlanarGeometryResult` ConformationResult holding actual ω angle,
  ring-flip state, sp2 pyramidalisation, ring-pucker phase.

The two pieces are paired as one atom-of-work concept ("characterise
each atom's chemistry+geometry context") but they sit on different
scopes per `feedback_object_model_scope_discipline`.

### F. No "named debt" carries forward

The IUPAC-topology revert's three-debts list (per-residue prochiral
verification, variant overrides, DFT-compare completeness) does NOT
carry into this work. The CIP verification debt dissolves into the
RDKit-primary methodology. The variant-overrides "debt" is just
part of the upcoming generalisation. DFT-compare completeness is
its own thread (memory entry
`project_dft_compare_calculator_completeness`), unrelated to the
substrate work.

### G. C++ generator, not Python (locked)

The generator binary is C++ linking cifpp + RDKit, not a Python
script using gemmi-Python + rdkit-Python. The user's reasoning:
"the class of errors we make in python extends to absolute loss of
information undetected far far more often." C++'s stricter typing
catches more silent bugs.

### H. Process log committed alongside the generated source

The user's instruction: "make sure we log process and inspect log,
even if it is put away once." The generation log
(`src/generated/LegacyAmberSemanticTables.log.txt`) is committed
alongside the typed-enum source (`.cpp`). The user inspects once
per regeneration; committed logs are the reproducibility audit
trail for future readers and downstream stats analysis.

---

## What's next (sequenced)

The remaining Session B work is genuinely separable from what
landed; it benefits from a fresh context. Order matters less than
the architecture-clarity that's already in place.

### Next-session priority 1: generalise synthesised-field tables to all 20 residues

`tools/topology/build_semantic_tables.cpp::SynthesisedForPhe` is the
template. The generalisation is:

1. **Add `SynthesisedFor<Residue>` per residue** (ALA, ARG, ASN,
   ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PRO, SER, THR,
   TRP, TYR, VAL). Each follows the PHE pattern:
   - Backbone atoms get `PlanarGroupKind::PeptideAmide`.
   - Ring atoms get the corresponding `RingSystemKind` and
     `RingPositionLabel` per the chemistry tables in the research
     dossier `spec/plan/topology-fields-research-2026-05-05.md`.
   - Polar Hs get the appropriate `PolarHKind`.
   - Pseudoatom-group members get the corresponding
     `PseudoatomMembership`.

2. **Variant handling** for HID / HIE / HIP / ASH / GLH / CYX / CYM /
   LYN / ARN / TYM. The CCD has separate entries for some (HID,
   HIE, HIP, ASH, GLH, CYM); for others (CYX, LYN, ARN, TYM) the
   variant is a derivation from the canonical (drop / add specific
   Hs, adjust formal charge). Use the CCD entry where it exists;
   derive otherwise. The variant-table from
   `src/AminoAcidType.cpp::AMINO_ACID_TYPES` documents the
   per-variant atom presence.

3. **Markley alternation for prochiral Hs.** RDKit's CIPLabeler
   labels chiral CENTRES (CA → S for L-amino acids), not prochiral
   Hs. The pro-R / pro-S mapping for diastereotopic Hs follows
   from Markley 1998 Figure 1: for L-amino acids, β-onward
   methylenes have Hβ3 = ProR, Hβ2 = ProS; Gly inverts (HA2 =
   ProR, HA3 = ProS). Encode as a small per-residue alternation
   rule that reads `DiastereotopicIndex` from the parser and emits
   `ProchiralStereo` per atom. The full per-residue table is in the
   research dossier "Field 3" section.

4. **Coverage test** (Task #10): every (residue, variant, atom) has
   all 14 fields populated with non-`None` provenance source. Lives
   at `tests/test_legacy_amber_semantic_tables.cpp`. Asserts
   programmatically that the generator's output is complete for the
   full standard 20.

### Next-session priority 2: integrate generated table into runtime build

The generated `src/generated/LegacyAmberSemanticTables.cpp` is
checked in but not yet compiled into `libnmr_shielding`. To
integrate:

- Add `src/generated/LegacyAmberSemanticTables.cpp` to the
  `nmr_shielding` source list in the main `CMakeLists.txt`.
- Optionally add a small lookup function in the same generated file
  or in a parallel `src/LegacyAmberSemanticTablesAccess.{h,cpp}`
  that, given `(AminoAcid, atom_local_idx, variant_idx)`, returns
  a `const AtomSemanticTable*` (or nullptr if out of range).

This is mechanical work; the design is already in
`spec/plan/topology-substrate-implementation-plan-2026-05-05.md`.

### Then: Session C+D as planned

After the generator is complete and the runtime can read its tables:

- Extend `LegacyAmberTopology` with the 14 typed fields populated
  via the lookup at `Protein::FinalizeConstruction`.
- New `PlanarGeometryResult` ConformationResult for the
  conformation-side ω, ring-flip, sp2 pyramidalisation,
  ring-pucker geometry per frame.
- H5 emission: `/atoms/legacy_amber/<14 fields>` plus
  `/conformation/planar_geometry/<angles>`.
- SDK catalog entries.
- Tests as per the implementation plan.

### Then: Session E (IUPAC + BMRB projections)

The Crystal Projection Rule consumers. Pure functions on
`LegacyAmberTopology`; output goes through H5 boundary as projected
name columns. User flagged E as "a must" for analysis to continue.

---

## File map

What lives where after this session.

### Spec / plan

```text
spec/plan/topology-fields-research-2026-05-05.md           [research dossier]
spec/plan/topology-substrate-implementation-plan-2026-05-05.md [the central spec]
spec/plan/session-handoff-20260505.md                       [this file]
```

### Runtime headers (consumed by both generator and the future populator)

```text
src/SemanticEnums.h     [the typed-enum vocabulary; 14 fields + provenance + AtomSemanticTable]
src/Types.h             [pre-existing; Hybridisation, AminoAcid, Element]
```

### Generator binary (build-time only; off by default)

```text
tools/topology/CMakeLists.txt          [opt-in subdir gated by NMR_BUILD_TABLE_GENERATOR=ON]
tools/topology/build_semantic_tables.cpp [the generator, ~1100 lines]
tools/topology/README.md                [how to regenerate]
```

### Generated artefacts (committed, not compiled into runtime yet)

```text
src/generated/LegacyAmberSemanticTables.cpp    [PHE only at end of this session]
src/generated/LegacyAmberSemanticTables.log.txt [structured generation log]
```

### Tests

```text
tests/test_string_barrier.cpp  [5 tests, all pass; backstops the four-layer barrier]
```

### Memory entries (auto-loaded next session)

```text
project_legacyambertopology_semantic_fields_a2.md   [the 14 fields are live forward intent]
feedback_planned_calculators_stay_planned.md        [calculator notes are not cruft]
project_iupac_topology_landed_20260426.md           [STALE BY REVERT marker]
feedback_aimnet2_required_no_weasel.md              [pre-existing, applicable]
```

---

## Pitfalls / footguns

The user has burned multiple times on these patterns; flagged here
so future-me does not repeat.

1. **Do NOT investigate the reverted IUPAC-topology episode.** The
   `iupac-fix-attempt-archive-2026-04-27.tar.gz` exists and must
   not be extracted. CLAUDE.md is explicit. The current substrate
   work is the architectural successor; the prior attempt's code
   is not relevant.

2. **Do NOT add RDKit or gemmi to runtime src/.** This includes:
   - Don't `#include "rdkit/..."` in any src/ header or .cpp.
   - Don't add an `nmr-runtime-rdkit` link target.
   - The string-barrier test will fail loudly if you do.
   - The discipline is enforced by the linker; bypassing the
     test-level enforcement would still fail at link time.

3. **Do NOT remove the process log.** The user explicitly wants the
   structured generation log committed alongside each regeneration.
   It is the audit trail; downstream stats analysis joins on it.
   "Put away once" means inspected once per regeneration, not
   deleted.

4. **Do NOT pre-frame the synthesised fields as "weaker provenance."**
   Per Decision D above. The chemistry citations make these fields
   first-class contributions.

5. **Do NOT split the substrate work from its conformation companion
   beyond the substrate-vs-conformation scope split.** The user
   pushed back on excessive subdivision: "you never met a
   distinction that doesn't justify another seperate phase or
   session, but sometimes you are right." The substrate fields and
   per-frame planar-geometry CR land paired in C+D, not on separate
   tracks.

6. **Do NOT change the C++ vs Python generator decision.** Locked
   per Decision G. C++ wins on type safety; Python's silent failure
   modes are the issue.

7. **Avoid `Hybridisation` typo.** The existing project enum is
   `Hybridisation` (British spelling) in `Types.h:91`. The
   generator + SemanticEnums.h follow this. Don't introduce
   `Hybridization` (US spelling) -- it would shadow.

8. **`row.get<T>(key)` in cifpp returns `std::tuple<T>`,** not
   `std::optional<T>`. Use `row[key].value_or<T>(default)` for
   single-value extraction with default. Pattern is documented
   inline in the generator.

9. **`numAromaticRings` does not exist on RDKit 2023.09's
   `RingInfo`.** Count manually by walking `bondRings()` and
   checking `getIsAromatic()` on each bond. Pattern is in
   `LogRdkitPerception`.

10. **CIPLabeler does NOT label prochiral Hs.** It labels chiral
    centres (CA → S for L-amino acids). Pro-R / pro-S for
    diastereotopic methylene Hs comes from Markley alternation
    rule applied AFTER CIPLabeler establishes the parent's
    chirality. This is the next-session work.

---

## How to verify state in a fresh session

Quick checks a fresh session can run to confirm the state described
above:

```bash
# 1. HEAD matches.
git log --oneline -6
# Expect to see 721e681, b3e281b, 6d4aa12, d915fb4, 7390a6d, beea8cf.

# 2. String-barrier tests pass.
cmake --build build --target string_barrier_tests
./build/string_barrier_tests
# Expect 5/5 pass.

# 3. Linker-level barrier intact.
nm build/libnmr_shielding.a | grep -c '_ZN5RDKit'
# Expect 0.

# 4. Generator builds and runs.
cmake -B build-gen -DNMR_BUILD_TABLE_GENERATOR=ON
cmake --build build-gen --target build_semantic_tables
./build-gen/tools/topology/build_semantic_tables \
    --ccd data/ccd/components.cif \
    --output /tmp/sm_test.cpp \
    --log /tmp/sm_test.log
# Expect "OK -- generated /tmp/sm_test.cpp with 1 residue table(s)."

# 5. Generated file is well-formed C++.
g++ -std=c++17 -I src -I /usr/include/eigen3 -c \
    src/generated/LegacyAmberSemanticTables.cpp -o /tmp/synch.o
# Expect zero warnings/errors.

# 6. Process log structure.
head -20 src/generated/LegacyAmberSemanticTables.log.txt
# Expect "==== run-info ====" / "==== ccd-load ====" / etc.
```

---

## Closing note

Five commits this session, all green. The architecture is proven on
PHE end-to-end. The string barrier is layered and enforced. The
process log is the audit trail. No named debt carried forward.

Next session generalises the synthesised tables to the remaining 19
residues + 10 protonation variants and adds the coverage test. That
is mechanical work against an established template -- the chemistry
per residue is the substantive content; the architecture is settled.

Then runtime integration (Session C+D), then projections (Session E).
The user identified E as "a must" for analysis to continue, so it is
not optional even if the substrate is fully landed.

The pacing concern the user flagged is real: each session has been
making substantive design progress that is hard to reconstruct from
git diffs alone. This handoff plus the implementation plan plus the
research dossier give a fresh session enough to pick up without
losing direction. If anything in this document seems unclear at
session-pickup time, that is itself a signal -- the right next move
is to read all three docs in this order: research dossier, then
implementation plan, then this handoff.
