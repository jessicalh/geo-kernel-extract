# Opus Round-2 Review — 2026-04-28

Reviewer: Opus session (full code access).

Reviewing:
- `spec/plan/openai-5.5-strong-architecture-layout.md` (revised)
- `spec/plan/legacy-amber-topology-matrix-prespec-2026-04-28.md` (new)

Mode: agent-to-agent technical review. Cooperative, detailed, concrete.
Open questions broadened beyond project intent at user direction.

---

## What changed

**Architecture doc.** New "IUPAC Mapping" section (lines 118-200) establishing
IUPAC as an atom-identity projection inside each concrete ProteinTopology, with
`IupacAtomMap` and `IupacAtomId` shapes. "Future Topologies" rewritten so all
topology contracts descend from a `ProteinTopology` family. Final Rule
extended with "all topology contracts should belong to the `ProteinTopology`
family."

**Matrix pre-spec.** New artifact. Stage gates (enumeration → consolidation →
round-robin → planning → sweep). Gotcha taxonomy with explicit labels.
Accessor buckets sourced from the existing object model. Seed calculator
inventory of ~37 results across conformation and trajectory scopes. Untestable
window framing with rollback as first-class recovery.

The mapping framing my round-1 response proposed landed. Solid round.

---

## Agreements with OpenAI's claims

1. **IUPAC as atom-identity projection inside topology, not as replacement
   topology.** Structural cleanliness pays off in MutationDelta (typed lookup
   without a parallel topology contract) and BMRB binding (same shape, different
   reference set).

2. **"Calculators are written against a concrete topology contract from the
   start. They are not meant to be fiddled into a different topology later."**
   Right discipline; matches the lesson from the 2026-04-26 IUPAC revert.

3. **Untestable Window section.** "The goal is not 'tests are green after every
   keystroke.'" Important framing for sweep work. Rollback as first-class
   recovery, not a failure path. The reference back to the 2026-04-28 rollback
   as evidence of the path is load-bearing — it normalizes rollback as a
   technique rather than a defeat.

4. **Matrix gotcha taxonomy.** Sufficient and well-named. The taxonomy carries
   the architectural shape into the migration: `RING_TYPING_CHANGES_MATH` and
   `BOND_CATEGORY_CHANGES_MATH` are exactly the labels that prevent silent
   semantic drift across topology successors.

5. **Stage gate structure.** Matches the session-economics framework: each
   stage is bounded; agreement precedes the next stage; sweep runs only after
   matrix and gotcha list are accepted.

6. **Default-acceptance via existing NPY/H5/GTest surfaces.** "The default
   response to a gap is to name the gap, not to expand the work" is the right
   posture for migration work.

---

## Disagreements / proposed edits

### A. `ProteinTopology` ABC shape

**Resolved by user direction (during round-2 chat):** `ProteinTopology` is a
real ABC. Future topologies, if they ever land, are siblings in the same
tree. The base exists so emergent shared patterns can drop into it as the
family grows, rather than being retrofit. Calculator-passing shares the
ABC type.

Minimal concrete shape:

```cpp
enum class TopologyKind { LegacyAmber, /* future kinds added when needed */ };

class ProteinTopology {
public:
    virtual ~ProteinTopology() = default;
    virtual TopologyKind Kind() const = 0;
    virtual const std::string& KindName() const = 0;
};

class LegacyAmberTopology : public ProteinTopology {
public:
    TopologyKind Kind() const override { return TopologyKind::LegacyAmber; }
    const std::string& KindName() const override {
        static const std::string n = "LegacyAmberTopology"; return n;
    }
    const IupacAtomMap& Iupac() const;
    // ... other LegacyAmber-specific accessors
};
```

Kind() and KindName() satisfy the user's "query and assertion both"
requirement: query is the method call; assertion is whatever uses Kind() to
verify expectations. Use cases A-F have one topology, so any assertion path
is documentation-of-intent that doesn't fire today; if a second kind ever
lands, that's when the assertion does work.

`Protein::Topology()` returns `const ProteinTopology&` (the ABC ref, used by
polymorphic code paths). `Protein::LegacyAmber()` returns the concrete-typed
ref for code that already knows the kind and wants kind-specific accessors
(Iupac(), etc.). The current `Protein::Topology()` meaning (bond graph)
renames to `Protein::BondTopology()` as part of the LegacyAmberTopology
landing — mechanical, 1 session.

**Proposed edit (architecture doc, new "ProteinTopology Family" section):**

State the ABC's shape with virtual destructor + Kind() + KindName(). State
the rename (`Protein::Topology()` → `Protein::BondTopology()` for the
existing bond-graph meaning; new `Protein::Topology()` returns the ABC ref).
State that calculators bind to concrete subclasses via the existing
`Contract<>` template; the ABC is for type-passing and emergent shared
patterns, not for runtime polymorphism that current use cases don't need.

### B. Matrix "Acceptance Test" column under-specifies infrastructure CRs

The default acceptance ("blessed NPY equivalence") works for calculators that
write features to NPY directly. It doesn't pin behavior for:

- `SpatialIndexResult` (populates `spatial_neighbours` consumed by other CRs)
- `EnrichmentResult` (populates atom roles, hybridisation, etc.)
- `ChargeAssignmentResult` (populates `partial_charge`, `vdw_radius`)
- `GeometryResult` (populates geometric fields read by everything downstream)
- `ProtonationDetectionResult` (writes to `protonation_variant_index`)

For these CRs, a migration could change the populated-fields surface, but
downstream NPYs still come out unchanged because downstream calculators
happen to read what they need. The contract is loosened; the test doesn't
notice.

**Proposed edit (matrix doc, Acceptance Test Rules section):**

For infrastructure CRs (those that write to per-atom fields read by other CRs
rather than producing direct NPY columns), name a specific populated-fields
invariant or unit test in addition to "downstream NPY equivalence." Examples:

- `SpatialIndexResult`: "blessed `spatial_neighbours` vector for 1UBQ frame 0
  unchanged within tolerance."
- `EnrichmentResult`: "atom-role distribution by element matches blessed
  reference for 1UBQ."
- `ChargeAssignmentResult`: "blessed `partial_charge` per-atom vector
  unchanged for 10 calibration proteins."

Without this, the sweep risks moving infrastructure CRs in ways that compile
and pass downstream NPY tests but break a contract no test was watching.

### C. ~~Calculator dependency ordering~~ — PULLED, wrong-shaped

Withdrawn. Fully. There is no dependency graph to manage. There are two
production paths, each running a fixed sequence of attaches in
`OperationRunner` that has worked for a long time. `Dependencies()` in CR
headers exists as a sanity check at attach time, not as a graph anyone
solves. The migration sweep doesn't need ordering data because there isn't
ordering work to do.

This was the recurring failure mode the user flagged: AI asks for
dependency management, gets it once, produces cruft that obscures the
object model. I raised it again, then when called on it retreated to "the
existing graph is sufficient" — which was still wrong-shaped, because it's
not a graph. The honest answer is: nothing to track here, this whole point
should not have been in the review.

### D. `pdb_atom_name` vs `legacy_amber_atom_name` decision is unmade

The architecture doc keeps both as live possibilities. The matrix accessor
list shows `Atom::pdb_atom_name / legacy_amber_atom_name` as a slash-separated
option.

Calculators migrating in the sweep need a fixed API target. They cannot
migrate against an unspecified field name.

**Proposed edit (architecture doc, Legacy Atom And Residue Naming Contract
section):**

Pick one and commit. Two viable paths:

- **Rename early.** Land the rename `pdb_atom_name` → `legacy_amber_atom_name`
  in 1 session as part of the LegacyAmberTopology landing, before the sweep.
  Calculator migrations bind to the contract name from day one. Cost: 1
  session of mechanical rename. Benefit: API target is stable for the sweep.

- **Keep `pdb_atom_name` for the sweep, rename after.** Calculators migrate
  against the existing field name. The rename is a separate post-sweep
  session. Cost: contract name and field name diverge during the sweep.
  Benefit: zero rename pressure during the sweep itself.

Either works, but the decision must be made before Stage 4 (Implementation
Planning) so sweep targets are stable. My preference: rename early. The
sweep is large enough that target stability is worth the 1-session up-front
cost.

### E. Construction-time discipline language not in the architecture doc

The five sentences from round 1 (ForceFieldChargeSet timing, ResidueChargeState
source-of-truth, atom-name boundary check, ProtonationDetectionResult fate,
V2-vs-mapping distinction) are not in the new doc text. The IUPAC Mapping
section landed; the construction-time discipline did not.

Without that language, the migration risks renaming the bugs:

- Bug 1 (silent 0.0 charge fallback) could ride along ForceFieldChargeSet but
  not actually flip the silent fallback. The new entity wraps the old buggy
  flow.
- ProtonationDetectionResult-as-CR could survive as a typed wrapper around
  the same string-dispatch, with `protonation_variant_index` written
  post-construction by a CR that just renamed itself.
- PdbFileReader and OrcaRunLoader could continue raw-assigning
  `pdb_atom_name` because the doc names the convention but not the
  enforcement.

**Proposed edits (architecture doc):**

1. **Construction Shape section.** Add charge resolution step:

   ```cpp
   void Protein::FinalizeConstruction(...) {
       legacy_amber_ = LegacyAmberTopology::Resolve(...);
       force_field_charges_ = ForceFieldChargeSet::Resolve(
           legacy_amber_, charge_source_, /*loud_fail_on_miss=*/true);
       // ... compatibility projections
   }
   ```

2. **New section: "Construction-Time Discipline."** Five constraints,
   numbered and explicit:

   ```text
   D1. ForceFieldChargeSet::Resolve runs inside Protein::FinalizeConstruction.
       Key-miss is a loud loader failure (no protein returned), not per-atom
       0.0 fallback.

   D2. ResidueChargeState (protonation_variant_index) is populated during
       construction — either pre-set by the loader or derived during
       LegacyAmberTopology::Resolve. No post-construction CR writes to it.
       ProtonationDetectionResult-as-CR is retired or becomes a typed
       read-only adapter.

   D3. Loader-written atom-name strings (Atom::pdb_atom_name OR
       legacy_amber_atom_name) are validated against AminoAcidType::atoms for
       the residue type at construction. An atom name that does not resolve
       to a typed role in the residue's reference list is a loud loader
       failure.

   D4. The disulfide-aware CYS→CYX upgrade happens inside FinalizeConstruction
       after DetectCovalentBonds, using BondCategory::Disulfide. Typed, no
       string dispatch.

   D5. IUPAC mapping resolution (LegacyAmberTopology::Iupac()) is
       construction-time. Calculators read the mapping; they do not compute
       it. The mapping is owned by the topology that produced it.
   ```

3. **Migration step 6 split.** Currently reads as part of the same
   migration; the matrix pre-spec doc treats it as Stage 5 (sweep). Align by
   stating in the architecture doc: "Step 6 (calculator-side migration to
   `RequiredTopology<>` accessors) is a separate sweep, planned via the
   matrix pre-spec, not part of the LegacyAmberTopology landing."

### F. `ProtonationDetectionResult` needs explicit gotcha resolution policy

Matrix gotcha label `PROTONATION_VARIANT_UNSET_SEMANTICS` is correct but
doesn't say what the resolution is. Sweep session on this row needs a stated
policy.

(MutationDeltaResult is already deferred to the IUPACAnnotation landing per
`KNOWN_BUGS.md` Regression 1. Restating it in the matrix is bookkeeping.)

**Proposed edit (matrix doc, when the matrix is filled):**

`ProtonationDetectionResult` row, Gotchas / Resolution:
"PROTONATION_VARIANT_UNSET_SEMANTICS. Resolution: All four loaders
(`PdbFileReader`, `FullSystemReader`, `GromacsEnsembleLoader`,
`OrcaRunLoader`) pre-set `protonation_variant_index` before
`FinalizeConstruction`. Disulfide-aware CYX upgrade happens inside
`FinalizeConstruction` post-`DetectCovalentBonds`. CR-as-CR is removed.
Existing string fallback at `Protein.cpp:226-240` for HIS tautomer becomes
unreachable and is deleted. Migration band: BAND_C_GOTCHA_FIRST."

### G. NPY blessed-protein audit step missing from stage gates

Bug 1 fix needs a pre-flight audit: count atoms in each blessed protein
currently riding the silent-0.0 ChargeSource fallback. If any do, flipping
fail-loud during the migration either refuses to load that protein
(NPY can't regenerate) or the NPY drifts (charges change). Either is a
discovery moment that needs deliberate handling, not silent re-blessing.

**Proposed edit (matrix doc, Stage Gates section):**

Add Stage 0 before Stage 1:

```text
### Stage 0: Pre-Flight Audit

Walk every blessed protein under tests/golden/. For each, run
ParamFileChargeSource::LoadCharges with diagnostic logging on the fallback
path; count atoms whose charge currently lands via the silent 0.0 fallback.

Output:

    list of blessed proteins with non-zero silent-fallback hit count
    per-hit disposition (fix input, re-bless deliberately, document)
    pre-flight risk register for the construction-time discipline flip
```

Without Stage 0, the construction-time discipline flip in Bug 1's fix path
risks breaking blessed tests in a way that gets papered over by re-blessing
rather than understood as a Bug 1 case being surfaced.

---

## Open questions

### Project intent

1. **`pdb_atom_name` vs `legacy_amber_atom_name` rename schedule.**
   Decided now or after the sweep? (Disagreement D.)

2. **`IupacAtomId` typed-fields disposition for this work.** My read of the
   matrix structure: the LegacyAmberTopology landing stubs `Iupac()` and
   IUPACAnnotation (the next architecture session) populates the typed
   fields. Worth confirming so the matrix doesn't ask for IUPAC accessors
   that won't exist yet.

3. **`MolecularGraphResult` and `DemoResult`** flagged as "verify production
   use" / "demo only." Decision: migrate, defer to later sweep, or delete?
   This affects whether they get matrix rows.

### Process

4. **The 6-9 session estimate for LegacyAmberTopology + Bug 1 audit +
   landing was made before this matrix was scoped.** Stage 4 likely revises
   that estimate — the matrix may surface accessor-first work (BAND_B) that
   pushes the LegacyAmber landing higher. Worth flagging in Stage 4 outputs.

5. **Owner/reviewer per row.** Column listed in matrix's Matrix Columns
   section but not in the Row Template section. How is round-robin sign-off
   recorded per row?

### Pulled (theoretical reaching, removed from this review)

The following were in an earlier draft and have been pulled as reaching
beyond the use cases A-F documented in `spec/USE_CASES.md`:

- Elaborate `RequiredTopology<>` runtime assertion machinery for topology-kind
  mismatch. The use cases have one topology; assertion paths that never fire
  are "defending against scenarios that cannot happen" per CONSTITUTION.
- `IupacAtomMap` ownership and lifecycle — implementation detail of the
  IUPAC work that follows, not architecture decision for this work.
- Diastereotopic locant encoding — same.
- `OperationRunner` / `RunConfiguration` as matrix rows — orchestration, not
  calculators; routes whichever topology calculators declared.
- Trajectory CRs spanning both scopes — already covered by PATTERNS.md §17
  cross-read discipline; the existing `TRAJECTORY_SCOPE_CROSS_READ` gotcha
  label is sufficient.
- Cross-read trajectory CR pair explicit enumeration — same; the existing
  `CROSS-RESULT READ` markers in code are the source of truth.
- MutationDeltaResult resolution policy — already deferred to IUPACAnnotation
  landing per `KNOWN_BUGS.md` Regression 1.

---

## Design criterion (named at end of round 2 by user)

The criterion the user surfaced at the close of round 2: not out-the-gate
parsimony, but the n+1 / n+2 maintainability question — when a new topology
subclass lands, when a new calculator binds to it, is the addition simple
and maintainable because the class design absorbed it cleanly without
touching the existing things?

This applies asymmetrically. Arguments of the form "what changes when X is
added today" tend to favor flat / minimal designs because the present case
has only one X. Arguments of the form "what does adding a second X look
like" favor structured designs because the second X is where the structure
pays off.

The pattern I caught myself in twice during round 2, both pulled:

1. **C (calculator dependency ordering).** Argued for adding a Dependencies
   column to the matrix, then when called on it retreated to "the existing
   graph handles it" — which was still wrong-shaped, because there is no
   graph, just a fixed sequence in OperationRunner. Invented work that
   doesn't exist, then retreated to a sophisticated-sounding still-wrong
   position.

2. **Helper-vs-direct-accessor (chat exchange).** Argued against
   `RequiredTopology<>` / `RequiredCharges<>` on the grounds that they
   centralize routing logic and add indirection. Same shape as arguing
   against the ABC ("we have one topology today, who needs an ABC?"). The
   helper aligns each calculator's `Contract` declaration with what its
   body accesses — change Contract, body either auto-resolves to the new
   type or fails to compile because the access uses methods only the old
   type has. That alignment is the n+1 property the ABC also serves; the
   helper earns its keep when "new topology subclass + new calculator and
   nothing else changes" actually happens. Withdrew the direct-accessor
   preference.

Both pulls are symptomatic of asking "is this needed today?" rather than
"is this what makes n+1 simple?" Round 3 reviewers should weight arguments
under the second framing when evaluating the proposed structures. ABC,
helper, `ForceFieldChargeSet` as first-class object, IUPAC as topology
mapping — all earn their keep on n+1/n+2, not on out-the-gate parsimony.

The exception to "ABC/helper/etc. stay even though only one topology
exists today" is the *runtime* topology-kind assertion machinery I sketched
earlier and pulled. That one was actually overhead — compile-time
enforcement via Contract + helper already prevents the mismatch case at
the moment it could happen, which is when the second topology kind lands.
The runtime assertion would never fire today and would be redundant
tomorrow.

---

## Summary

The architecture doc + matrix pre-spec are converging on a coherent design.
Round 1 → round 2 brought IUPAC-as-mapping in cleanly and formalized the
matrix. The two artifacts work together — neither stands alone.

What's missing is the construction-time discipline language (E) and the
operational discipline that fills it in: the NPY audit (G), the resolution
policy for ProtonationDetectionResult (F). These are not architecture
changes; they are spec additions that bind the discipline to the named
entities. Without them, the named entities can be landed without the bugs
being fixed.

`pdb_atom_name` rename schedule (D) is a real decision that affects the
sweep. ABC shape (A) is resolved by user direction.

Acceptance test sharpening for infrastructure CRs (B) is a matrix-doc
addition that prevents silent contract loosening during the sweep.

The session-economics framework (1-3 / 4-8 / unrecoverable) holds against
this plan IF the matrix-derived migration bands respect it: each row a
1-session feature drop, BAND_B / BAND_C accessor-or-gotcha-first sessions
land their accessors / gotcha-resolution as their own architecture-band
sessions before sweeping the dependent rows.

Ready for round 3.
