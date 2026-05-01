# Opus Round-3 Review — 2026-04-28

Reviewer: Opus session (full code access).

Reviewing:
- `spec/plan/openai-5.5-strong-architecture-layout.md` (round-3 revision)
- `spec/plan/legacy-amber-topology-matrix-prespec-2026-04-28.md` (round-3 revision)
- `spec/plan/README.md` (round-3 revision)

Mode: agent-to-agent technical review.

---

## What changed in round 3

**Architecture doc.**

- New "Required ABC Object Model" section. ABC declared explicitly with
  `Kind() / Name() / AtomCount() / ResidueCount() / IupacOrNull()` virtuals.
  `LegacyAmberTopology` shown as `final` concrete subclass. Explicit rule:
  "Do not use runtime polymorphism to make calculators topology-agnostic.
  There is no planned switch-up where an old calculator is handed a richer
  topology and expected to keep meaning the same thing."
- New "Result And Table Naming" section. `*Result` is reserved for
  ConformationResult / TrajectoryResult. Charge entities use `*ChargeTable`.
  `ForceFieldChargeTable`, `MopacChargeTable`, `EeqChargeTable`,
  `Aimnet2ChargeTable`, `CalculatedChargeTable`, `NoChargeTable`.
- New Protein API split: `TopologyBase()` returning the ABC ref,
  `TopologyAs<TopologyT>()` returning the concrete-typed ref,
  `LegacyAmber()` retained as compatibility shorthand.
- Helper signature now enforces ABC membership at compile time:
  `static_assert(std::is_base_of_v<ProteinTopology, TopologyT>)`.
- New ownership names on Protein: `protein_topology_`,
  `force_field_charges_`, `legacy_bond_topology_projection_`. Explicit
  rule: do not keep using a member named `topology_` for the covalent
  graph.
- Construction Shape updated: `protein_topology_ = LegacyAmberTopology::Resolve(...)`.
- Explicit statement on contract-axis orthogonality: "topology contract
  and charge-table contract are intentionally orthogonal. That
  orthogonality matters because this project is expected to need at least
  one more concrete topology before the maths work is done."

**Matrix doc.**

- New "ProteinTopology ABC" accessor bucket, deliberately small.
- Matrix Columns now include `Required topology concrete class` and
  `Needed ProteinTopology ABC accessors`.
- Stage 3 output adds `accepted ProteinTopology ABC shape`.
- Stage 4 output adds `ProteinTopology ABC landing order`.

**README.**

- New "Current baseline" section: prior half-finished reviewer edits are
  not active instructions; respond to the packet.
- New "Review protocol" section with five rules, including: "Do not
  replace the ABC topology model with an untyped mapping model."

---

## Agreements (round-2 concerns that landed)

1. **ABC made explicit and concrete.** `Kind() / Name() / AtomCount() /
   ResidueCount() / IupacOrNull()` is the right minimal surface. The
   `final` keyword on `LegacyAmberTopology` and the prohibition on
   topology-agnostic calculators forecloses the "calculators portable
   across topologies" failure mode that would have undone the n+1
   property.

2. **Calculators bind concrete, not ABC.** The doc states this in three
   places (Required ABC, Calculator Contract Convention, Final Rule). The
   `static_assert` in `RequiredTopology<>` enforces it at compile time. No
   slippage path for a future calculator to accidentally bind to the base.

3. **Result/Table naming distinction.** Real catch — `Set` was overloaded
   with set-of-atoms in normal usage; `Table` is concrete (rows keyed by
   atom_index). Resolves the naming clash with the `*Result` suffix.

4. **Orthogonality made explicit.** "Two template axes are intentionally
   orthogonal" is the right framing for the n+1 case. Explicit statement
   "this project is expected to need at least one more concrete topology"
   matches the user's design criterion.

5. **`Protein::TopologyBase()` + `TopologyAs<TopologyT>()` split.** Names
   the polymorphic and concrete entry points distinctly. New code uses
   the concrete-typed accessor; polymorphic ref available when needed.

6. **Anti-drift README guardrail.** "Do not replace the ABC topology
   model with an untyped mapping model" — pinning it as a review-protocol
   rule means subsequent reviewers can't quietly back out of the ABC
   commitment.

7. **`IupacOrNull()` nullability — resolved by user.** It's a transition
   state: at LegacyAmberTopology landing, no concrete topology has an
   IUPAC mapping yet; the ABC accessor returns null until IUPACAnnotation
   populates it. Same transition shape works for any future topology
   landing before its own IUPAC mapping.

---

## Outstanding from round 2 (still concrete)

The structural model is solid. What hasn't moved across rounds is the
operational discipline that determines whether the migration fixes the
existing bugs or just renames them. These are not architecture-doc
decisions — they're spec sentences that need to live somewhere.

### O1. Construction-time charge resolution policy

Architecture doc's Construction Shape (lines 866-887) shows topology
resolution but no charge resolution. Migration step 1 reads:

> Add `ForceFieldChargeTable` and make `ChargeAssignmentResult` expose it
> while continuing to populate `ConformationAtom::partial_charge` and
> `ConformationAtom::vdw_radius`.

This naturally reads as wrap-the-existing-flow. If `ForceFieldChargeTable`
is built lazily by the CR rather than at construction, and the
lookup-miss policy isn't specified, Bug 1's silent-0.0 fallback survives
the migration.

The grounded data: `ff14sb_params.dat` covers all 20 standard residues
plus variants HIE/HID/HIP/ASH/GLH/CYX/LYN. Standard fallback rows for
HIS/ASP/GLU catch variant-extra hydrogens via secondary lookup. Real miss
conditions for the four loaders working correctly are: protein has CYM /
ARN / TYM (rare deprotonated states), or input PDB had non-Standard atom
names (user-error, fixable by translation at boundary). Under standard
production paths the miss rate is essentially zero. Per user direction:
"if it is very unlikely then we fail."

**Proposed home:** matrix row for `ChargeAssignmentResult`. Specifically
the row's gotcha-resolution column (which doesn't exist yet — see E1
below).

**Proposed sentence:** "ForceFieldChargeTable is built during
Protein::FinalizeConstruction (or in the loader, before the loader
returns). On `(residue, legacy_amber_atom_name)` lookup miss against
ff14sb_params.dat, the loader refuses to build the protein with an error
naming the offending pair. No per-atom 0.0 silent fallback. NamingRegistry
retry is available if a real workflow surfaces the convention-mismatch
case."

### O2. ProtonationDetectionResult retirement policy

Architecture doc names `protonation_variant_index` as part of
LegacyAmberTopology's residue symbolic topology and says it's
"load-bearing," but doesn't say who populates it, when, or what happens
to the existing `ProtonationDetectionResult`-as-CR. Today the CR runs
post-construction and does the seven post-construction string-dispatch
sites identified in `spec/EVIL_STRING_AUDIT_2026-04-28.md`.

**Proposed home:** matrix row for `ProtonationDetectionResult`, gotcha
`PROTONATION_VARIANT_UNSET_SEMANTICS`, resolution column.

**Proposed sentence:** "All four loaders (PdbFileReader,
FullSystemReader, GromacsEnsembleLoader, OrcaRunLoader) pre-set
protonation_variant_index before FinalizeConstruction. Disulfide-aware
CYX upgrade happens inside FinalizeConstruction post-DetectCovalentBonds
(typed via BondCategory::Disulfide, no string dispatch). The CR-as-CR is
removed. The string fallback for HIS tautomer at Protein.cpp:226-240
becomes unreachable and is deleted. Migration band:
BAND_C_GOTCHA_FIRST."

### O3. NPY blessed-protein Stage 0 audit

Bug 1's fail-loud discipline (O1) needs a pre-flight: count atoms in
each blessed protein currently riding the silent-0.0 fallback. If any
are found, surface as Bug 1 cases to address explicitly (fix input,
re-bless deliberately, document) before the discipline flips. Without
this, the migration risks breaking blessed tests in a way that gets
papered over by re-blessing rather than understood as a Bug 1 case being
surfaced.

**Proposed home:** matrix Stage Gates section, new Stage 0 before
Stage 1.

**Proposed text:**

```text
### Stage 0: Pre-Flight Audit

Walk every blessed protein under tests/golden/. For each, run the
existing ParamFileChargeSource::LoadCharges with diagnostic logging on
the fallback path; count atoms whose charge currently lands via the
silent 0.0 fallback. Same shape for any other construction-time discipline
that is about to flip from silent to loud.

Output:

    list of blessed proteins with non-zero silent-fallback hit count
    per-hit disposition (fix input, re-bless deliberately, document)
    pre-flight risk register for the construction-time discipline flip
```

### O4. `pdb_atom_name` vs `legacy_amber_atom_name` decision

Architecture doc still keeps both as live possibilities. Matrix accessor
list still shows `Atom::pdb_atom_name / legacy_amber_atom_name`. This is
an API-target decision that needs to be made before the calculator sweep
begins, otherwise migrating calculators bind to an unspecified field
name.

Two viable paths:
- **Rename early** as a 1-session mechanical change before the sweep
  starts. Sweep targets are stable.
- **Keep current name through the sweep, rename after.** Field name and
  contract name diverge during the sweep.

**Proposed home:** Stage 4 (Implementation Planning) decision. Surface in
Stage 3 round-robin agreement.

---

## New, smaller items from this round

### N1. Matrix row template lacks a Resolution field per gotcha

The gotcha taxonomy is exact (POST_CONSTRUCTION_STRING_IDENTITY,
CHARGE_SOURCE_MATCHING, etc.) but the row template has only "Gotchas"
as a field. A label without a resolution policy is half a fact —
someone migrating that row in a sweep session needs to know what the
fix is, not just that one is needed.

**Proposed edit (matrix doc, Row Template):**

Add a `Gotcha resolutions` field below `Gotchas`. Format: one line per
gotcha label, naming the resolution policy. This is where O1, O2, and
the future-bundled MutationDelta resolution land per-row.

### N2. `TopologyAs<TopologyT>()` semantics not shown

The ABC-helper section shows `TopologyAs<TopologyT>()` being called from
inside `RequiredTopology<>` but doesn't specify the cast semantics. With
one concrete topology this doesn't matter; once a second concrete topology
lands, "what does TopologyAs do when the actual topology is a different
concrete kind than requested?" becomes a real call.

Three options:
- Static cast (assumes caller knows the kind; UB if wrong)
- Dynamic cast (runtime check, returns null/throws on mismatch)
- Kind-guarded static cast (check `Kind() == expected`, then static_cast)

The Kind() field on the ABC is suggestive of the third. Worth pinning
in the architecture doc when a second topology lands; not a decision
that affects this work.

### N3. Owner/reviewer column gap in matrix Row Template

Mentioned in round 2; still present. Column listed at line 158, not in
the Row Template at lines 451-481. Mechanical fix.

---

## Open questions

### Project intent

1. **Where do O1 / O2 / O3 sentences land — architecture doc or matrix
   row resolutions?** My proposal above puts them in the matrix because
   they are per-calculator behavior commitments rather than universal
   type-level rules. But if "fail loud at construction on charge miss"
   is meant as a system-wide rule rather than a per-row policy, it
   belongs in the architecture doc. The user's earlier "the rest lands"
   suggested these are agreed in principle; the open question is where
   they live as text.

2. **`pdb_atom_name` rename schedule.** Decided in this work or after?

### Process

3. **Round-robin agreement gate for the architecture doc.** The doc has
   gone through two reviewer-driven revisions (rounds 2 and 3). At what
   point does it close — after a fourth round with no objections? After
   explicit user sign-off? Some signal that "the doc is now baseline,
   matrix-fill begins" is needed before Stage 1 (Enumeration) starts.

---

## Net read

Architecture is solid after round 3. ABC + concrete subclass + orthogonal
contract axes + ChargeTable family + IUPAC as topology mapping with
nullable transition + small ABC surface. The n+1 case (a richer concrete
topology landing as a sibling subclass while existing calculators stay
on LegacyAmberTopology) is supported by the structure as written.

What needs to happen before implementation begins is the operational
discipline. Three sentences (O1, O2, O3) — one resolving Bug 1's
fail-loud policy, one retiring ProtonationDetectionResult as a CR, one
pre-flight audit before flipping discipline. Plus the rename decision
(O4). Plus the matrix template gap (N1) so per-row resolutions have a
home.

Then the matrix gets filled across the seed inventory, gotchas get their
resolutions named, untestable window gets defined, and Stage 5 sweep
begins as bounded 1-session calculator migrations against blessed NPY
acceptance.

Ready for round 4.
