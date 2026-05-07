# Slice A pre-commit audit (2026-05-07, adversarial)

Read-only audit of the substrate-side Slice A changes — the typed
surface for `ProPyrrolidineRing` + Pro per-atom labels + TRP indole
9-atom perimeter (tertiary slot) + the aromatic-ring-type-count named
constant. Per the contract: Slice A is substrate-only, no calculator
code touched, no iteration-surface change, no NPY ABI shift.

Scope of audit: every diffed file, generated table consistency,
predicate rename safety, test coverage of the new typed surface,
sacred-doc consistency, h5-reader / SDK ordinal-compat documentation,
exhaustive-switch behaviour. Pre-existing `SmokeTest.WithDft` not
flagged per audit instructions.

---

## Executive summary

Slice A is **mostly clean and ready for commit**, with **0 CRITICAL**
findings, **6 CONCERN** findings, and **5 NITPICK** findings.

The substrate changes are structurally sound: the additive
`tertiary` slot is consistently emitted across all 503 generated rows
(1509 RingMembership values total, exactly 3 per row), Pro 5-label
substrate emits correctly through the cap-merge path (the existing
`ProNTermPreservesPyrrolidineRing` test confirms preservation under
NTERM cap-delta), the new `ProPyrrolidineRing` class exhausts its
seven required virtuals, the `CreateRing` factory is now exhaustive
with explicit `Count: break;` + abort (good fail-loud discipline),
and the predicate rename (`InAnyRing` → `IsInAnyRing`) has no
external callers — the rename is fully contained.

The CONCERN findings are mostly about coverage gaps and subtle
documentation drift the user typically catches:

1. **The `static_assert` on `kAromaticRingTypeCount == 8` is
   tautological** — it pins the constant against its own defining
   literal, not against any external invariant. Anyone who edits
   the literal `8` to `9` simultaneously updates the assert. The
   comment claims the assert "pins the constant to that ABI value
   so that any future change to the constant fails the build until
   the ABI / NPY migration is also done"; in fact the assert provides
   no such pin.
2. **No test exercises the tertiary-slot encoding for TRP.** The
   substrate now flags 15 atoms (9 heavy + 6 H) with
   `Indole_Trp_9 / PerimeterMember`, but there is no test that
   reads `LegacyAmber().SemanticAt(ai).ring_position.tertiary.ring`
   for a TRP fixture. The encoding is structurally consistent but
   semantically untested.
3. **No test exercises the new substrate predicates**
   (`HasPrimaryRing`, `HasSecondaryRing`, `HasTertiaryRing`,
   `MembershipCount`, `IsPopulated`, `IsInAnyRing`). They are
   `constexpr` so the compiler verifies they compile; no test
   verifies their semantic correctness against fixture data —
   particularly that `MembershipCount` returns 3 for TRP
   bridgeheads CE2/CD2 and 2 for non-bridgehead perimeter atoms.
4. **`OBJECT_MODEL.md` quotes the old enum** (8 values, `Count = 8`)
   at lines 1054-1064. Slice A modifies the enum and OBJECT_MODEL.md
   is sacred-doc per `feedback_sacred_docs` memory entry, but isn't
   in the diff.
5. **h5-reader / Python SDK have stale "ordinal-compatible"
   documentation comments** that reference now-shifted line numbers
   in `src/Types.h` and don't surface that the new
   `RingTypeIndex::ProPyrrolidine = 8` and `RingAromaticity::None`
   are not represented downstream. Library will never serialise
   these to H5 / NPY today (calculator-side adoption is deferred),
   so this is documentation drift, not a runtime issue. But the
   "ordinal-compatible" claim is now strictly false.
6. **The TRP residue table in the residue reference doc does NOT
   include a `ring_tertiary` column** even though every perimeter
   atom has a populated tertiary slot. The schema bullet says the
   column is omitted only when "NotInRing for every atom in the
   residue" — but for TRP this condition is false (every perimeter
   atom has Indole_Trp_9 in tertiary). The prose bullet after the
   table describes the encoding correctly; the table itself is
   incomplete.

The NITPICK findings are wording / consistency cleanups that don't
affect correctness.

**Verdict: clean to commit; the static-assert finding (#1) is the
strongest CONCERN — easy to fix in this slice (one line) and
mismatches the documented intent loud-and-clear. The rest are
cleanly deferrable.**

---

## Findings table

| Severity | File:line | One-line description |
|---|---|---|
| CONCERN | `src/Types.h:227-231` | `static_assert(kAromaticRingTypeCount == 8, ...)` is tautological; pins the constant to its own defining literal, not to ABI / NPY shape |
| CONCERN | `tests/topology/test_legacy_amber_semantic_integration.cpp` (no test) | No test exercises `RingPosition.tertiary` slot population on TRP fixtures |
| CONCERN | `tests/test_ring_hierarchy.cpp` + tests/topology (no test) | New predicates `HasPrimaryRing`/`HasSecondaryRing`/`HasTertiaryRing`/`MembershipCount`/`IsPopulated`/`IsInAnyRing` have zero test coverage |
| CONCERN | `OBJECT_MODEL.md:1054-1064` | Sacred doc quotes old `RingTypeIndex` enum (8 values, Count=8) — Slice A diff doesn't update it |
| CONCERN | `h5-reader/src/model/Types.h:236, 240-259` + `python/nmr_extract/_types.py:6-17` | Downstream "ordinal-compatible with nmr::RingTypeIndex" claims are stale; no flag-down on the now-asymmetric enum |
| CONCERN | `spec/plan/topology-residue-reference-2026-05-05.md:715-733` | TRP residue table does not include `ring_tertiary` column even though every perimeter atom has populated tertiary slot |
| NITPICK | `src/Types.h:233-246` | `RingTypeName` switch retains `default: return "?";` — masks future "missing case" warnings; inconsistent with `Ring.cpp::CreateRing`'s explicit `Count: break;` discipline |
| NITPICK | `tools/topology/build_semantic_tables.cpp:1706` | Comment "minus the shared bridgehead bond" is technically the bond, not the bridgehead atoms; bridgeheads CE2/CD2 ARE in the perimeter |
| NITPICK | `tools/topology/build_semantic_tables.cpp:2515-2528, 2530-2553` | `RingSystemLiteral` and `RingPositionLabelLiteral` switch statements have unreachable trailing `return ::nmr::RingSystemKind::NotInRing;` — exhaustive switches over enum classes don't need them; can mask future missing-case warnings |
| NITPICK | `spec/meta-docs-review/CONSTITUTION_CRITIQUE.md:374` | Stale claim "All 8 ring types as classes" — now 9 |
| NITPICK | `tests/topology/test_legacy_amber_semantic_integration.cpp:594-636` | `MakeNTermProProtein` fixture builds Pro at NTERM with H1/H2/H3, but Pro NTERM_CHARGED has only 2 NH protons biologically (secondary amine becomes +1 ammonium with 2 Hs). Outside Slice A's scope but worth flagging since the test is the primary Pro-cap-merge regression gate |

---

## Per-finding deep dives

### CONCERN-1: `static_assert(kAromaticRingTypeCount == 8, ...)` is tautological

**File:** `/shared/2026Thesis/nmr-shielding/src/Types.h:226-231`

**Evidence:**

```cpp
// Line 226: declaration
inline constexpr int kAromaticRingTypeCount = 8;
// Lines 227-231: static_assert
static_assert(kAromaticRingTypeCount == 8,
              "kAromaticRingTypeCount is the per-aromatic-type NPY ABI "
              "shape; changing it requires migrating the ConformationAtom "
              "per_type_* arrays + the calculator emission code + the "
              "downstream NPY schema together.");
```

**Why it matters:**

The comment immediately above this assert (lines 218-225) says:

> The static_assert below pins the constant to that ABI value so that
> any future change to the constant fails the build until the ABI /
> NPY migration is also done.

But `static_assert(X == 8)` where `X` is defined on the line above as
`= 8` is checking the constant against its own literal. Anyone who
changes line 226 from `= 8` to `= 9` will simultaneously update line
227 from `== 8` to `== 9` (or the assert will fail in a way that's
trivially fixable by editing one line, with no flag to anything else).

The intent is clear — pin the constant to the load-bearing literal
`8` that appears in `ConformationAtom.h:131-134, 244-249` and
`fileformat/analysis_file.h:99-102, 162-170` — but no compile-time
relation exists between the named constant and the literals at those
sites.

A real pin would be one (or more) of:
- `static_assert(static_cast<int>(RingTypeIndex::ProPyrrolidine) == kAromaticRingTypeCount, "...")` —
  this links the constant to the boundary it names (Pro is the first
  saturated ring); changing one without the other is caught.
- A `static_assert` IN `ConformationAtom.h` or `fileformat/analysis_file.h`
  that the literal `8` equals `kAromaticRingTypeCount` — this links
  the constant to the actual ABI sites.

The current assert (a) gives the *appearance* of pinning while not
actually pinning, (b) makes the comment misleading.

**Recommended action:** Add a meaningful pin to ConformationAtom.h
or replace the tautological assert in Types.h with one against
`static_cast<int>(RingTypeIndex::ProPyrrolidine)`. (The audit doesn't
prescribe which.)

---

### CONCERN-2: No test exercises the TRP `tertiary` slot population

**File:** `/shared/2026Thesis/nmr-shielding/tests/topology/test_legacy_amber_semantic_integration.cpp`
(no test added in Slice A diff)

**Evidence:**

Grep for tertiary slot uses across the entire repo (`tests/`,
`src/`, `h5-reader/`, `ui/`, `python/`):

```
src/SemanticEnums.h:683:    RingMembership tertiary;
src/SemanticEnums.h:691, 693, 700: predicate definitions
tools/topology/build_semantic_tables.cpp:1710-1714: generator
tools/topology/build_semantic_tables.cpp:2599: emission
src/generated/LegacyAmberSemanticTables.cpp: 488 default + 15 populated
src/generated/LegacyAmberSemanticTables.h:120-131: comment
spec/plan/topology-encoding-dependencies-2026-05-05.md: doc
spec/plan/topology-residue-reference-2026-05-05.md: doc
```

No test file references `ring_position.tertiary` or `Indole_Trp_9`.
The 15-atom TRP perimeter substrate emission is structurally
consistent (verified: 1509 = 503 * 3 RingMemberships) but no test
asserts `SemanticAt(CG_atom).ring_position.tertiary.ring ==
Indole_Trp_9` for a TRP fixture, nor that all 9 heavy + 6 H atoms
share the `Indole_Trp_9 / PerimeterMember / 9 / arom=true /
n_heteroatoms=1` value.

**Why it matters:**

The substrate is the source of truth for Bundle B's
`ConstructRingsFromSubstrate`. If Bundle B reads the tertiary slot
to detect TRP perimeter membership, a regression in the generator
that emits `NotInRing` for tertiary on TRP atoms would be caught
only at Bundle B integration, not at Slice A regen. Bundle B is
likely 2+ commits later. A simple substrate-dump assertion in Slice
A would close this gap.

The new typed surface (RingPosition.tertiary, Indole_Trp_9,
PerimeterMember) is the surface Bundle B depends on. Slice A leaves
this surface untested.

**Recommended action:** flag for test coverage; engineering lead
decides whether to add a TRP perimeter test in Slice A or defer to
Slice B's first consumer.

---

### CONCERN-3: New predicates `HasPrimaryRing`/`HasSecondaryRing`/`HasTertiaryRing`/`MembershipCount`/`IsPopulated`/`IsInAnyRing` have zero test coverage

**File:** `/shared/2026Thesis/nmr-shielding/src/SemanticEnums.h:691-705,
926-928` (definitions); no test file exercises them.

**Evidence:**

```
$ grep -rn "MembershipCount\|HasPrimaryRing\|HasSecondaryRing\|HasTertiaryRing" \
       --include="*.cpp" --include="*.h" 2>/dev/null

src/SemanticEnums.h:691-705: definitions only
```

No callers anywhere in the codebase outside the definitions
themselves. The constexpr nature means they compile correctly, but
their *semantic* correctness is unverified:

- Does `MembershipCount` return 3 for TRP CE2/CD2 (bridgeheads,
  primary=Indole_Trp_5, secondary=Indole_Trp_6, tertiary=Indole_Trp_9)?
- Does `MembershipCount` return 2 for TRP non-bridgehead perimeter
  atoms (e.g. CG: primary=Indole_Trp_5, secondary=NotInRing,
  tertiary=Indole_Trp_9)?
- Does `MembershipCount` return 1 for non-TRP ring atoms (e.g. PHE
  CG: primary=Benzene_Phe, others NotInRing)?
- Does `MembershipCount` return 0 for backbone atoms?
- Does `IsPopulated` correctly distinguish `NotInRing` (default)
  from any populated slot?

The substrate tables emit data that *should* satisfy these, but the
predicates themselves have no test that verifies they read the
substrate correctly.

The audit-doc-comment at SemanticEnums.h:686-689 explicitly notes
the non-bridgehead-perimeter case ("a TRP non-bridgehead perimeter
atom carries {primary=Indole_Trp_5 (or Indole_Trp_6),
secondary=NotInRing, tertiary=Indole_Trp_9} — primary AND tertiary
populated, with secondary empty. Use these predicates to ask the
right question, not count-based predicates.") — but no test verifies
this is what the substrate actually emits.

**Recommended action:** flag for test coverage gap; engineering lead
decides.

---

### CONCERN-4: `OBJECT_MODEL.md:1054-1064` quotes old enum

**File:** `/shared/2026Thesis/nmr-shielding/OBJECT_MODEL.md:1051-1065`

**Evidence:**

```
1054	enum class RingTypeIndex {
1055	    PheBenzene    = 0,
1056	    TyrPhenol     = 1,
1057	    TrpBenzene    = 2,
1058	    TrpPyrrole    = 3,
1059	    TrpPerimeter  = 4,
1060	    HisImidazole  = 5,
1061	    HidImidazole  = 6,
1062	    HieImidazole  = 7,
1063	    Count         = 8
1064	};
```

This is the OLD (8-value) enum. The Slice A change updates
`src/Types.h` to add `ProPyrrolidine = 8, Count = 9`. OBJECT_MODEL.md
is not in the diff and remains stale.

**Why it matters:**

Per `feedback_sacred_docs` memory entry: "GEOMETRIC_KERNEL_CATALOGUE,
OBJECT_MODEL, CONSTITUTION: read every session, no exceptions." The
sacred docs are session-startup reading. A future session reading
OBJECT_MODEL.md will see 8 values and be misled.

Additionally, the new substrate names (`ProPyrrolidineRing`,
`RingAromaticity::None`, `RingSystemKind::Indole_Trp_9`,
`RingPositionLabel::PerimeterMember`, the five Pro labels, and the
extended `RingPosition.tertiary` slot) are NOT in OBJECT_MODEL.md or
PATTERNS.md anywhere:

```
$ grep -n "ProPyrrolidine\|ProRing\|kAromaticRingTypeCount\|Indole_Trp_9\|tertiary" \
       OBJECT_MODEL.md PATTERNS.md spec/CONSTITUTION.md
(no output)
```

**Recommended action:** Update OBJECT_MODEL.md (at minimum the
RingTypeIndex enum quote at line 1054 plus a paragraph on the
9-value enum / Pro Ring / 3-slot RingPosition); engineering lead
decides whether sacred-doc updates land in Slice A's commit or as a
follow-on.

---

### CONCERN-5: Stale "ordinal-compatible" docs in h5-reader and Python SDK

**Files:**
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/model/Types.h:230-259`
- `/shared/2026Thesis/nmr-shielding/python/nmr_extract/_types.py:1-17`

**Evidence (h5-reader):**

```cpp
// h5-reader/src/model/Types.h:233
// Ordinal-compatible with nmr::RingAromaticity (src/Types.h:168).
enum class RingAromaticity { Full, Reduced, Weak };

// h5-reader/src/model/Types.h:240-259
// RingTypeIndex — 8 concrete ring types.
// Ordinal-compatible with nmr::RingTypeIndex (src/Types.h:175..185).
enum class RingTypeIndex : int32_t {
    PheBenzene    = 0,
    ...
    HieImidazole  = 7,
};
constexpr int RingTypeCount = 8;
```

After Slice A:
- `nmr::RingAromaticity` has 4 values (`Full, Reduced, Weak, None`),
  not 3. h5-reader has 3.
- `nmr::RingTypeIndex` has 9 values, h5-reader has 8.
- `src/Types.h:175..185` line range no longer matches; the enum is
  now at `src/Types.h:175..199`.

**Evidence (Python SDK):**

```python
# python/nmr_extract/_types.py:1
"""Physical enums matching the C++ RingTypeIndex and BondCategory."""

class RingType(IntEnum):
    PHE = 0
    ...
    HIE = 7

N_RING_TYPES = 8
```

The Python SDK's `RingType` IntEnum has 8 values; library has 9.
The docstring "matching the C++ RingTypeIndex" is now strictly
false (8 vs 9). The constant `N_RING_TYPES = 8` matches
`kAromaticRingTypeCount`, not `RingTypeIndex::Count` — but the docstring
doesn't distinguish.

**Why it matters:**

This is intentional under the locked Path D scope: calculator-side
adoption is deferred, so the library never serialises ring_type=8
(ProPyrrolidine) to NPY/H5, and downstream readers never see the
new value. **The runtime is safe.** But:

- Documentation drift: comments now lie. Line numbers in
  h5-reader's Types.h ("src/Types.h:175..185") are wrong. The
  "ordinal-compatible" claim is partial (0-7 still match; 8 does
  not exist downstream). New eyes reading these comments will be
  misled.
- The asymmetry is part of the locked design but is not
  surfaced as such anywhere downstream. A Bundle B / per-calculator
  follow-on slice that does adopt ring_type=8 in NPY emission would
  break h5-reader's `DecodeRingType` (rejects n=8 since `n < 8`
  fails) and the Python SDK's `RingType(8)` (raises `ValueError`).
  The h5-reader's loader has explicit handling
  (`h5-reader/src/io/QtProteinLoader.cpp:99-107`) that returns
  `RingTypeIndex::PheBenzene` as a sentinel on out-of-range, with
  an error message — graceful, but misleading: PHE is real
  chemistry, not a sentinel.

**Recommended action:** Optional — Slice A could include a comment
in h5-reader's Types.h and python/nmr_extract/_types.py noting the
deferred-extension convention ("library has 9 values; downstream
only sees 0-7 because Pro Ring is not yet emitted"). But this is
genuinely deferred work and adding a comment now risks committing
to a specific extension shape. Engineering lead decides.

---

### CONCERN-6: TRP residue table omits `ring_tertiary` column

**File:** `/shared/2026Thesis/nmr-shielding/spec/plan/topology-residue-reference-2026-05-05.md:715-733`

**Evidence:**

The schema bullet at lines 101-110 says:

> - `ring_primary` / `ring_secondary` / `ring_tertiary`:
>   `RingMembership` shorthand `sys/pos/size/arom/het` where: ...
>   - The `ring_tertiary` slot is omitted from per-residue tables
>     when it is `NotInRing` for every atom in the residue (the
>     common case). It is populated for the indole 9-atom perimeter
>     on TRP (`RingSystemKind::Indole_Trp_9` per Case 1995, J.
>     Biomol. NMR 6, 341-346) — the conjugated π current circuit
>     encoded as a typed substrate slot alongside the chemical
>     5-ring + 6-ring decomposition.

But the TRP residue table at lines 715-733 has only 2 columns
(`ring_primary | ring_secondary`):

```
| atom | element | planar_group | planar_stereo | pseudoatom | polarH | ring_primary | ring_secondary | prochiral | formal_chg | exch | notes |
|---|---|---|---|---|---|---|---|---|---|---|---|
| N | N | PeptideAmide | NA | — | NotPolar | NotInRing | NotInRing | ...
...
| CG | C | Aromatic5Ring | NA | — | NotPolar | Indole_Trp_5/Ipso/5/t/1 | NotInRing | ...
| CD1 | C | Aromatic5Ring | NA | — | NotPolar | Indole_Trp_5/PyrroleBeta/5/t/1 | NotInRing | ...
...
```

For TRP, the schema bullet's condition for omitting (`NotInRing for
every atom`) is FALSE — every perimeter atom has a populated
tertiary slot. The table should either include the `ring_tertiary`
column or the schema bullet should be reworded to allow prose-only
description.

The notes section at line 729 has a prose paragraph describing the
perimeter encoding correctly — but readers who just look at the table
miss the data. The doc is inconsistent with itself.

**Why it matters:**

The residue reference doc is the canonical encoding spec the
generator implements against. A future regen pass that fixes a
generator bug by re-reading the doc could miss the perimeter
encoding because the doc's TRP table doesn't show it.

**Recommended action:** Either add `ring_tertiary` column to the TRP
table OR reword schema bullet 109 to acknowledge prose-only
description. Engineering lead decides.

---

### NITPICK-1: `RingTypeName` switch retains `default: return "?";`

**File:** `/shared/2026Thesis/nmr-shielding/src/Types.h:233-246`

**Evidence:**

```cpp
inline const char* RingTypeName(RingTypeIndex t) {
    switch (t) {
        case RingTypeIndex::PheBenzene:     return "PHE";
        ...
        case RingTypeIndex::ProPyrrolidine: return "PRO";
        default: return "?";
    }
}
```

Compare with the new (better) `Ring.cpp::CreateRing` discipline
(lines 52-75 of Ring.cpp): explicit `case Count: break;` for the
sentinel, no `default:` clause, fall-through to abort with
diagnostic. Adding a new `RingTypeIndex` value to that switch will
fail compilation with `-Werror=switch` (or warn with `-Wswitch`).

`RingTypeName`'s `default:` swallows new enum values silently,
returning `"?"` instead of forcing an explicit case. Inconsistent
with `CreateRing`'s discipline.

**Recommended action:** Optional cleanup — replace `default:` with
explicit `case RingTypeIndex::Count: break;` and a final `return
"?";` after the switch (or abort). Not load-bearing.

---

### NITPICK-2: Generator comment "minus the shared bridgehead bond"

**File:** `/shared/2026Thesis/nmr-shielding/tools/topology/build_semantic_tables.cpp:1702-1707`

**Evidence:**

```
1702    // Indole 9-atom perimeter (Case 1995, J. Biomol. NMR 6, 341-346):
1703    // a third ring system covering all 9 heavy atoms of the conjugated
1704    // π current circuit. Encoded in the tertiary RingMembership slot.
1705    // The perimeter is chemistry-derived from the union of the 5-ring
1706    // and 6-ring atoms minus the shared bridgehead bond (CE2-CD2).
1707    // Ring-attached H atoms inherit perimeter membership per the
1708    // convention used for primary/secondary slots.
```

"The shared bridgehead bond (CE2-CD2)" is the bond between the two
bridgehead atoms — but the perimeter walk INCLUDES both CE2 and
CD2 (they're at the junction of the fused system; perimeter
traversal turns the corner at each). What's "removed" from the
perimeter is the CE2-CD2 bond (which is the shared edge of the two
rings), not the bridgehead atoms themselves.

The 5-ring has atoms `{CG, CD1, NE1, CE2, CD2}` (5 atoms; closing
edge CD2-CG). The 6-ring has atoms `{CD2, CE2, CZ2, CH2, CZ3, CE3}`
(6 atoms; closing edge CE3-CD2 wait — the 6-ring's closing bond is
CD2-CE3, and the bridge-bond is CE2-CD2). Union: 5+6 = 11; minus 2
shared atoms (CE2 + CD2 counted once): 9 unique atoms. Perimeter
walk: `CG → CD1 → NE1 → CE2 → CZ2 → CH2 → CZ3 → CE3 → CD2 → CG`
(closes back) — 9 atoms, traverses the OUTER boundary, never crosses
the CE2-CD2 bridge bond.

**Why it matters:**

Trivial wording. The substrate tagging is correct (all 9 heavy +
6 H atoms in the if-statement at lines 1717-1722). Only the
explanatory comment is sloppy.

**Recommended action:** Optional — reword to "minus the shared edge
bond between the bridgeheads" or similar. Not load-bearing.

---

### NITPICK-3: Trailing returns after exhaustive switch in generator helpers

**File:** `/shared/2026Thesis/nmr-shielding/tools/topology/build_semantic_tables.cpp:2515-2528, 2530-2553`

**Evidence:**

```cpp
const char* RingSystemLiteral(RingSystemKind r) {
    switch (r) {
        case RingSystemKind::NotInRing:       return "...";
        ... // all 8 enum values
        case RingSystemKind::Indole_Trp_9:    return "...";
    }
    return "nmr::RingSystemKind::NotInRing";  // line 2526
}
```

The switch covers all enum values explicitly (no `default:`). The
trailing return after the closing brace is unreachable. This is the
common GCC-pleasing idiom, but it ALSO means `-Wswitch` won't fire
if a new enum value is added without adding a case (the function
"returns" via the trailing return). This is the same masking issue
as NITPICK-1 above, applied to two different switches.

Same situation in `RingPositionLabelLiteral` at lines 2530-2553.

**Recommended action:** Optional — replace trailing `return ...;`
with `__builtin_unreachable()` or a fail-loud diagnostic. Or simply
keep current style; -Wswitch usually warns inside switch on
enum-class even without `default:` in the configurations gcc/clang
ship.

---

### NITPICK-4: `spec/meta-docs-review/CONSTITUTION_CRITIQUE.md:374` says "All 8 ring types"

**File:** `/shared/2026Thesis/nmr-shielding/spec/meta-docs-review/CONSTITUTION_CRITIQUE.md:374`

**Evidence:**

```
374    **TRUE.** All 8 ring types as classes with virtual properties.
```

Now stale (9 ring types). This is a meta-review doc, not a primary
spec, so consequence is low. But it's textual reality drift.

**Recommended action:** Optional — update to "9 ring types" or "8
aromatic ring types + 1 saturated".

---

### NITPICK-5: `MakeNTermProProtein` fixture uses 3 NTERM Hs on Pro

**File:** `/shared/2026Thesis/nmr-shielding/tests/topology/test_legacy_amber_semantic_integration.cpp:594-636`

**Evidence:**

```cpp
const std::pair<const char*, nmr::Element> atom_set[] = {
    {"N",   nmr::Element::N},
    ...
    {"H1",  nmr::Element::H},
    {"H2",  nmr::Element::H},
    {"H3",  nmr::Element::H},   // <-- third NTERM H on PRO
};
```

Pro at NTERM_CHARGED biologically has 2 NH protons (secondary amine
becomes +1 ammonium with N-CD bond, N-CA bond, and 2 N-H bonds —
total 4 substituents on tetrahedral N). H3 should not exist on Pro
NTERM_CHARGED.

The cap table for `kCapNtermCharged` at the generated lines 655-660
has 4 atoms (N, H1, H2, H3) — designed for NON-PRO residues. PRO at
NTERM should use a different cap variant (e.g. only N + H2 + H3 or
similar).

**Why it matters:**

The test passes today because:
1. Slice A is preserving existing behaviour, not introducing this fixture.
2. The substrate-test asserts only `ring_position.primary.{ring,
   position}` and `formal_charge` — not whether H3 has biologically
   valid Pro chemistry.

The test fixture pre-existed Slice A; the fixture's biological
correctness is outside Slice A's scope. **Just flagging because**
the test is the regression gate for the "PRO at NTERM preserves
ProRingNitrogen" property after Slice A's label swap, and an
unsound fixture undermines confidence in the regression gate.

**Recommended action:** Outside Slice A's scope. Possibly worth a
follow-on cap-table cleanup, but not blocking.

---

## What looked right

Many things in Slice A are well-executed and worth calling out
explicitly to calibrate the negative findings:

1. **`CreateRing` factory is now exhaustive with explicit `Count:
   break;` and `std::abort()` on fall-through.** The previous silent
   fallthrough to `PheBenzeneRing` was a real "blur" issue. The
   replacement is fail-loud, exhaustive (`-Wswitch` will fire if a
   new enum value is added without an explicit case), and includes
   a stderr diagnostic with the bad cast value. Excellent
   discipline.

2. **The predicate rename from `InAnyRing` / `InTwoRings` /
   `InThreeRings` to `HasPrimaryRing` / `HasSecondaryRing` /
   `HasTertiaryRing` / `MembershipCount` / `IsInAnyRing` is fully
   contained.** Grep confirms zero callers of the old names anywhere
   (production code, tests, viewer, h5-reader, Python SDK). The
   rename is safe. The new names are also semantically more honest
   (count-based predicates were broken for non-bridgehead perimeter
   atoms; slot-based predicates are accurate).

3. **The structural extension of `RingPosition` from `{primary,
   secondary}` to `{primary, secondary, tertiary}` is consistently
   applied.** Verified: 503 atom rows × 3 RingMemberships = 1509
   total RingMembership instances in the generated cpp. No partial
   rows. C++ aggregate-initialization handles the trailing default
   correctly because all generated rows now emit explicit
   `{primary, secondary, tertiary}` triples (the generator was
   updated to emit the 3rd slot, not relying on default-init).

4. **Cap-merge invariant preserved.** The `ApplyCapDelta` comment
   block (lines 115-131 of generated header) now correctly
   references `Pyrrolidine_Pro/ProRingNitrogen` instead of the
   stale `Pyrrolidine_Pro/Saturated`. The cap-merge logic itself
   was unchanged in Slice A (correctly — it preserves
   `ring_position` from chain regardless of label). The
   `ProNTermPreservesPyrrolidineRing` test confirms the new
   `ProRingNitrogen` survives the cap overlay.

5. **`ProPyrrolidineRing` class implements all seven required
   virtuals.** `Intensity()`, `LiteratureIntensity()`,
   `JBLobeOffset()`, `NitrogenCount()`, `Aromaticity()`, `RingSize`
   (inherited from `FiveMemberedRing`), `TypeName()`. Literal `0.0`
   for ring-current-related virtuals with citation (Joule & Mills
   2010 ch. 7) — physics, not calibration. Good.

6. **`ProPyrrolidineProperties` test pins all the new properties.**
   `EXPECT_DOUBLE_EQ` on Intensity, LiteratureIntensity,
   JBLobeOffset; `EXPECT_EQ` on NitrogenCount, Aromaticity,
   RingSizeValue. Comprehensive.

7. **`AromaticTypeCountBoundary` test pins the conceptual link
   between `kAromaticRingTypeCount` and `RingTypeIndex::ProPyrrolidine`**
   (`EXPECT_EQ(static_cast<int>(RingTypeIndex::ProPyrrolidine),
   kAromaticRingTypeCount)`). This is the link the static_assert in
   Types.h *should* have but doesn't. The test catches the issue
   at runtime; the static_assert is supposed to catch it at
   compile-time but currently doesn't.

8. **The Pro NTERM cap-merge regression test was correctly updated
   to expect the new `ProRingNitrogen` label**, and the failure
   message is detailed and cites the relevant chemistry literature
   (Vega & Boyer 1979, Schubert et al. 2002).

9. **Generator is exhaustive on its own switches**:
   `RingSystemLiteral` covers all 8 RingSystemKind values
   (including new `Indole_Trp_9`); `RingPositionLabelLiteral` covers
   all 20 RingPositionLabel values (original 14 + 5 new Pro labels +
   PerimeterMember). No new enum value is missed in the generator.

10. **`RingTypeName` includes the new `PRO` mapping.** Slice A
    correctly added the `ProPyrrolidine: return "PRO"` case.

11. **Documentation diff in `topology-residue-reference-2026-05-05.md`
    Pro section is comprehensive.** New 5-label encoding documented,
    H-inheritance convention documented, Joule & Mills citation for
    `Intensity = 0`, Vega & Boyer + Schubert citations for puckering
    chemistry, Locant/RingPosition orthogonality reaffirmed for Pro
    Cα.

12. **Documentation diff in `topology-encoding-dependencies-2026-05-05.md`
    Section A.3** lists every Slice A enum addition (RingAromaticity::None,
    RingSystemKind::Indole_Trp_9, RingPositionLabel Pro values,
    RingPositionLabel::PerimeterMember) with citations.

13. **No regression in existing aromatic ring substrate.** Spot-checks
    of HIS variants, PHE, TYR, TRP 5-ring + 6-ring slot populations
    show the additive structural extension preserves all previous
    chemistry (just adds default `NotInRing` tertiary on every atom
    that didn't get a populated tertiary).

14. **All 118 ring-related and topology-integration tests pass**
    (`ctest -R "Ring|Pro|Trp|test_topology"` returns 100%, 118/118).

---

## What I did NOT verify

In the interest of full disclosure, here's scope I didn't fully cover:

1. **NPY ABI bit-identity with prior runs.** I did NOT regenerate
   the smoke-test fixtures and binary-diff against the bless. Per
   `project_smoke_test_bless_deferred_20260424` memory, the user
   has decided per-commit bless is not the gate; this audit accepted
   that. If Slice A's substrate change had perturbed
   `ConformationAtom`'s `per_type_*` arrays, the calculator code is
   unchanged so no perturbation should occur — but verification is
   "by argument" not "by binary diff".

2. **Building the project.** I did not run `cmake --build`; I trust
   the test results from `ctest` (which depends on the binary
   building successfully). The test summary confirms ctest 402/403
   pass.

3. **Per-residue diff comparison of the regenerated generated cpp.**
   The diff stat reports `1006 insertions / 590 deletions` in the
   1500+ line generated file. I spot-checked PRO and TRP rows but
   did NOT diff every other residue's tertiary-slot extension to
   confirm zero semantic regression on non-affected chemistry. The
   structural counts (1509 = 503 × 3) confirm uniformity, but a
   silent generator bug that, e.g., reordered RingMembership fields
   wouldn't be caught by counts alone.

4. **Trajectory / integration tests via real PDB / TPR fixtures.**
   I ran ring-specific tests (118/118 pass) but did not run the full
   suite (would have taken too long, and the prompt says ctest =
   402/403). I did not confirm this audit on a multi-protein
   integration fixture.

5. **Disagreement between this audit and the existing audit
   reports** in `spec/plan/ring-investigation-2026-05-06/`. I
   cross-read those reports for context but did not fully trace
   every claim they make against current HEAD; some statements
   there may be stale relative to what Slice A actually landed.

6. **Symbol table inspection of the built `libnmr_shielding.a`.**
   The prompt says `nm libnmr_shielding.a | grep _ZN5RDKit` returns
   0 — I trust this; did not re-run.

7. **Effects on ui/ subproject.** I checked h5-reader and Python
   SDK for ordinal-compat issues; did not check `ui/`. ui/ links
   the library directly; if it consumes `RingTypeIndex` it would
   see `Count = 9`, but I did not verify whether ui/ has any
   compile-time link or runtime read.

8. **GeometryChoiceBuilder records emitted by future Pro Ring
   iterations.** Slice A doesn't construct Pro rings (calculator
   side deferred); no GeometryChoice impact yet. Slice B's first
   consumer will exercise this path. Outside Slice A scope per the
   contract.
