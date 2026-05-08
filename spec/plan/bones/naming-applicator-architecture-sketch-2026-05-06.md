# Naming-application architecture — initial sketch (2026-05-06)

**Status:** WIP architecture for runtime atom-name canonicalisation.
This document is the agreed marching orders for the `NamingRegistry`
refactor that lands as the next bundle. It is NOT canonical patterns
yet — emergent issues found during implementation will reshape
specifics; PATTERNS.md and OBJECT_MODEL.md carry pointers here, not
generalised pattern descriptions, until the implementation lands and
Session E exercises the same shape across IUPAC + BMRB + AMBER
projections.

The architecture is described in algorithmic and concrete-data-model
terms. No religious adherence — common sense applies; consult when
unclear.

---

## Why this exists

The current `NamingRegistry` flat-table architecture has the
following manifest failures, surfaced by the 1Z9B fleet variance and
the ILE δ-methyl shift cases during Bundle B:

1. **Not idempotent.** Per-atom rules can rewrite a canonical input
   (e.g. PRO `HD2` with siblings `{HD2, HD3}`) to a non-canonical
   form (`HD3`, colliding with the existing `HD3`). Today this is
   masked by upstream consistency we don't enforce: the fleet
   inputs always have non-canonical names. A future load path that
   produces canonical names corrupts silently.

2. **Source-overloaded.** `ToolContext::Charmm` is the historic name
   of the source-context tag, but the rules tagged with it are
   AMBER `pdb2gmx`-RTP deviations, not CHARMM force-field naming.
   The enum lies about what the rules represent.

3. **No per-atom audit trail.** The current registry returns a single
   string output without recording which rule fired or why. A future
   triage agent asking "why did this atom canonicalise to X?" has to
   re-derive the answer from code inspection.

4. **Cross-rule-set conflicts hidden.** When two rules from different
   conceptual sources (e.g. an AMBER convention rule and a Pdb2gmx
   deviation rule) both target the same atom, the current
   architecture either silently picks one (first-match) or chains
   them (Stage 1 then Stage 2), without ever expressing which
   project decision authorises the choice.

5. **Self-documentation ad hoc.** Rule citations live as comments in
   `NamingRegistry.cpp`; there is no programmatic API to ask "what
   sources contribute to this canonicalisation?"

The `RecanonicaliseAfterProtonation` follow-up pass added in Bundle B
is a workaround that addresses (1) for the LYN HZ chain by detecting
duplicates after the chain rewrites — which is a guard that cancels
side effects of earlier rules. The sketched architecture replaces
that by making rules' predicates examine their context and only fire
when applicable.

---

## Algorithm

For each atom under canonicalisation:

**Step 1 — collect.** A containing object iterates every rule it
holds. For each rule, the rule's predicate `Applies(ctx)` is
evaluated against the atom's context. Every rule whose predicate
returns `true` contributes a record to a transient per-atom map.
Multiple rules from multiple sources may contribute records; the
map preserves all of them.

**Step 2 — resolve.** A method on the containing object takes the
per-atom map plus the context and returns the chosen output. The
method body is the explicit, documented choice point. Each branch
of the method names which rule sets it resolves and cites the
project decision authorising the choice. Possible outcomes:

- Map empty AND input matches a known canonical form for `ctx`:
  return input unchanged (idempotent on canonical inputs).
- Map empty AND input does not match any canonical form: fail-loud
  with a diagnostic naming the input, the residue, the source, the
  full context, and a note that no rule applies and the input is
  unrecognised.
- Map has exactly one record: return that record's proposed output.
- Map has multiple records: branch into the documented per-case
  resolution logic. Each branch unit-tested. If a combination fires
  that no branch covers: fail-loud (a new project decision is needed
  before this case can be resolved).

The map is **transient**. It is built per `Apply(ctx)` call,
consumed by the resolver, and discarded. The chosen output is what
persists (written onto `Atom.pdb_atom_name` by the caller). The
containing object MAY support a debug-logging mode that emits the
map plus the resolution decision to a side channel for triage; this
is opt-in, off by default, never on `Atom`.

**Step 3 — context and rules order are properties of the containing
object.** The containing object owns:

- A flat list of rules, each tagged with a `NamingSource` enum value.
- The resolution method body (with documented per-case branches).
- The canonicality oracle: given a `(residue, variant, terminal)`
  context, what set of atom names is canonical AMBER ff14SB form?
  Used by the resolver to detect "input already canonical" and skip
  rules that would shift it. The oracle reads from the substrate's
  per-(residue, variant) atom inventory.
- Optional logging-mode flag.

The containing object lives in its own module (initially: extends
`src/NamingRegistry.{h,cpp}`; rename to a more descriptive name
discussable at review). It is **not** on `Protein`. It is **not** on
`LegacyAmberTopology`. It is **not** on `Atom`. It is its own thing,
constructed once at process start (or per-`FinalizeConstruction`,
implementation choice), and called by loaders + by
`Protein::FinalizeConstruction`.

---

## Concrete object model (initial)

```cpp
namespace nmr {

// ============================================================================
// NamingSource — enum tag on each rule, identifying which scientific or
// project source authored the rule. Symmetric with SemanticSource on the
// substrate side. New sources = new enum values.
// ============================================================================

enum class NamingSource : uint8_t {
    Unknown                  = 0,
    AmberFf14SBCanonical     = 1,   ///< Project-internal canonical AMBER ff14SB
                                    ///< (the target of all canonicalisation).
    Pdb2gmxAmberRtpDeviation = 2,   ///< AMBER pdb2gmx writes side-chain methylenes
                                    ///< and methyl-bearing carbons in older AMBER
                                    ///< RTP convention; deviations from canonical.
    CifppPdbInput            = 3,   ///< cifpp's PDB parsing yields IUPAC convention
                                    ///< names (mostly canonical AMBER post-Markley).
    OrcaEcho                 = 4,   ///< ORCA NMR-output echoes the input atom names.
    CharmmLegacy             = 5,   ///< Quarantined CHARMM/XTC path (retired
                                    ///< 2026-05-02; rules retained for tests).
    Markley1998              = 6,   ///< Markley et al. 1998 J. Biomol. NMR 12:1-23
                                    ///< nomenclature recommendations.
    BmrbAtomNomTbl           = 7,   ///< BMRB nomenclature table (bmrb.io/ref_info).
    IupacIub1969             = 8,   ///< IUPAC-IUB 1969 tentative rules.
    ProjectSynthesis         = 9,   ///< Project-internal synthesised conventions
                                    ///< (e.g. Trp 6-ring perimeter labels per
                                    ///< topology-encoding-dependencies §C.3).
};

// ============================================================================
// NamingContext — the input record passed to Apply(). Carries everything a
// rule's predicate might examine. New context fields can be added as new rule
// types demand.
// ============================================================================

struct NamingContext {
    NamingSource source = NamingSource::Unknown;
    std::string  input_name;
    AminoAcid    residue_type      = AminoAcid::Unknown;
    int          variant_index     = -1;                       ///< -1 = unresolved
    TerminalState terminal_state   = TerminalState::Internal;

    // Sibling atom names IN THE INPUT FORM. Snapshot of the residue's atom
    // names at the start of the current pass through the applicator, before
    // any rules rewrite them. The snapshot is critical for shift-pair rules
    // (e.g. Pdb2gmx ILE HD1/HD2/HD3 ↔ canonical HD11/HD12/HD13): the rule's
    // predicate examines the snapshot to determine non-canonical state, NOT
    // the partially-renamed state.
    std::set<std::string_view> sibling_input_names;

    // Parent heavy atom's input name (for hydrogen disambiguation). Empty
    // for non-hydrogen atoms or when parent is unknown.
    std::string parent_input_name;

    // Diagnostics-only fields (not consumed by rule predicates; included in
    // fail-loud messages for triage).
    int         residue_sequence_number = 0;
    std::string chain_id;
};

// ============================================================================
// NamingRule — single rule. Has a source tag, an Applies predicate, an
// Output transform, plus name + rationale for diagnostics and self-docs.
// ============================================================================

struct NamingRule {
    NamingSource     source;
    std::string_view name;       ///< Stable identifier for diagnostics + tests.
    std::string_view rationale;  ///< One-line citation of the source decision.

    // Predicate: does this rule apply to this context?
    std::function<bool(const NamingContext&)> applies;

    // Output: given the context, what canonical form does this rule propose?
    std::function<std::string(const NamingContext&)> output;
};

// ============================================================================
// NamingApplication — one rule firing on one atom. Per-atom records of these
// are accumulated in the transient map.
// ============================================================================

struct NamingApplication {
    const NamingRule* rule;             ///< Which rule fired.
    std::string       proposed_output;  ///< What that rule proposes.
};

// ============================================================================
// NamingApplicator — the containing object. Owns rules, the resolver, the
// canonicality oracle. Lives in NamingRegistry.{h,cpp}; replaces the
// existing flat-map architecture.
// ============================================================================

class NamingApplicator {
public:
    NamingApplicator();   // populates rule list + sets up canonicality oracle

    // Single-atom canonicalisation. Builds the per-atom map, calls the
    // resolver, returns the canonical output.
    std::string Apply(const NamingContext& ctx) const;

    // Whole-residue canonicalisation: snapshots sibling names ONCE for the
    // residue, then iterates atoms. Necessary for shift-pair rules to read
    // the original sibling set, not the partially-rewritten set.
    //
    // input_names parallel to residue.atom_indices; output canonical names
    // parallel to input. Caller writes the outputs onto Atom.pdb_atom_name.
    std::vector<std::string> ApplyResidue(
        const std::vector<std::string>& input_names,
        AminoAcid residue_type,
        int variant_index,
        TerminalState terminal_state,
        NamingSource source,
        int residue_sequence_number,
        std::string_view chain_id) const;

    // Optional debug-logging mode (off by default).
    void SetDebugLogging(bool enabled);

private:
    // Step 1: iterate rules, return all applications.
    std::vector<NamingApplication> Collect(const NamingContext& ctx) const;

    // Step 2: resolve the per-atom map to a single output. Body is the
    // explicit per-case decision logic; each branch documented with a
    // citation to the project decision authorising the choice.
    std::string Resolve(const std::vector<NamingApplication>& applications,
                        const NamingContext& ctx) const;

    // Canonicality oracle: is `ctx.input_name` already a valid canonical
    // AMBER ff14SB atom name for `ctx.residue_type` + `ctx.variant_index`
    // + `ctx.terminal_state`? Reads the substrate's per-residue/per-variant
    // chain inventory. Used by Resolve to detect idempotent (already-canonical)
    // input.
    bool IsCanonical(const NamingContext& ctx) const;

    // Diagnostic emit for fail-loud paths. Includes ctx, the application
    // map, and a pointer to this sketch document.
    [[noreturn]] void FailUnresolved(
        const NamingContext& ctx,
        const std::vector<NamingApplication>& applications,
        std::string_view reason) const;

    std::vector<NamingRule> rules_;
    bool debug_logging_ = false;
};

// Process-wide singleton accessor. Constructed once.
NamingApplicator& GlobalNamingApplicator();

}  // namespace nmr
```

The shape choices above are starting points, not commitments:

- `std::function` per rule is one option; concrete subclasses of an ABC are
  another. The lighter-weight `std::function` form keeps rule definitions
  inline and readable; the heavier ABC form gives compile-time type safety.
  Either works. The agent picks; review at codex round.
- Singleton vs per-call construction: singleton is simpler; per-call gives
  per-protein customisation. Default singleton; raise as emergent issue if
  per-protein customisation surfaces.
- The canonicality oracle could be hand-coded (per-residue switch over atom
  names) or could read directly from the substrate (`AminoAcidType::atoms`
  + variant tables). Reading the substrate is cleaner; doing so requires
  the applicator to have a reference to the substrate (or the substrate to
  expose canonical-name-set queries). Open design question.

---

## Data flow at `Protein::FinalizeConstruction`

The applicator is invoked twice during construction, with progressively
richer context:

```text
Loader (PdbFileReader / FullSystemReader / OrcaRunLoader):
    For each atom in residue:
        ctx.source       = (loader's known source)
        ctx.input_name   = (raw atom name from source)
        ctx.residue_type = (resolved residue type)
        ctx.variant_index = -1                    // unresolved at load
        ctx.terminal_state = (computed from chain position)
        ctx.sibling_input_names = (snapshot of all atoms' raw names in residue)
        ctx.parent_input_name = (heavy-atom parent's raw name, if H atom)
        canonical = applicator.Apply(ctx)
        atom->pdb_atom_name = canonical

Protein::FinalizeConstruction (after CovalentTopology::Resolve and second-pass
ResolveProtonationStates):
    For each residue with protonation_state_resolved && variant_index >= 0:
        snapshot input_names from current atom.pdb_atom_name across the residue
        For each atom in residue:
            ctx.source       = (preserved from load — still recorded on residue)
            ctx.input_name   = atom.pdb_atom_name  // post-Pass-1 form
            ctx.variant_index = residue.protonation_variant_index
            ctx.terminal_state = residue.terminal_state
            ctx.sibling_input_names = snapshot
            // ...other fields unchanged...
            canonical = applicator.Apply(ctx)
            atom->pdb_atom_name = canonical
```

Two passes through the same applicator. Same `Apply(ctx)` call. Pass 1
rules whose predicates examine source-context fire; pass 2 rules whose
predicates require `variant_index >= 0` fire. The two-pass shape is
emergent from rule predicates examining the context, not architectural.

There is no "reprotonate mode" — the system runs through
`FinalizeConstruction` once per protein; the two passes are sequential
events in that single construction.

`RecanonicaliseAfterProtonation` (added in Bundle B) is replaced by
this. The applicator handles both passes uniformly.

---

## Worked example 1: 1Z9B LYS-labelled-LYN-chemistry

Fleet input `1Z9B/input.pdb` labels residue 28 as `LYS` and provides
side-chain hydrogens HZ1, HZ2 (no HZ3). Chemistry is LYN (neutral
amine NH2); naming is non-canonical (LYS naming conventions, not LYN).

### Pass 1 — at load

For atom HZ1 in residue 28 (LYS):

```text
ctx.source                  = NamingSource::Pdb2gmxAmberRtpDeviation
ctx.input_name              = "HZ1"
ctx.residue_type            = AminoAcid::LYS
ctx.variant_index           = -1                  // unresolved
ctx.terminal_state          = TerminalState::Internal
ctx.sibling_input_names     = { "N", "H", "CA", "HA", "CB", ..., "NZ", "HZ1", "HZ2" }
ctx.parent_input_name       = "NZ"
```

Rule iteration:

- Rule `LysAmmoniumHzPassThrough` (source `AmberFf14SBCanonical`):
  Applies if residue is LYS, input is HZ1/HZ2/HZ3, AND siblings contain all
  three (canonical charged LYS state). Predicate evaluates `false` — siblings
  don't contain HZ3.
- Rule `LynHz1ToHz2_NonCanonicalLynShift` (source `AmberFf14SBCanonical`):
  Applies if residue is LYS, variant_index unresolved or = LYN, input is HZ1,
  AND siblings contain HZ1+HZ2 but not HZ3 (signaling non-canonical LYN
  state). Predicate evaluates `true`. Output: "HZ2".

Per-atom map after collection: 1 record.

Resolver: 1 record → return its output. Apply returns "HZ2".

For atom HZ2 (same residue, same pass):

```text
ctx.input_name              = "HZ2"
ctx.sibling_input_names     = { ..., "HZ1", "HZ2" }   // SNAPSHOT — original names
```

- `LysAmmoniumHzPassThrough`: predicate `false` (no HZ3 in siblings).
- `LynHz2ToHz3_NonCanonicalLynShift`: predicate `true` (siblings indicate
  non-canonical LYN). Output: "HZ3".

Map: 1 record. Resolver returns "HZ3".

After Pass 1: residue 28 atoms have pdb_atom_name HZ2 (was HZ1) and HZ3
(was HZ2). The other LYS-relevant atoms (N, NZ, etc.) pass through unchanged
because no rule with a positive predicate applies, and they match canonical.

### Between passes

`Protein::ResolveProtonationStates(true)` runs. Sees residue 28 has HZ2
and HZ3 but no HZ1. Existing inference: variant_index = 0 (LYN).
`residue.protonation_state_resolved = true`.

### Pass 2 — post-protonation

For atom that was HZ1 (now HZ2 after Pass 1) in residue 28 (now known
LYN):

```text
ctx.source              = NamingSource::Pdb2gmxAmberRtpDeviation   // preserved
ctx.input_name          = "HZ2"      // current pdb_atom_name
ctx.residue_type        = AminoAcid::LYS
ctx.variant_index       = 0          // LYN — resolved
ctx.terminal_state      = TerminalState::Internal
ctx.sibling_input_names = { ..., "HZ2", "HZ3" }   // snapshot of post-Pass-1
```

Rule iteration:

- Rule `LynCanonicalHzPassThrough` (source `AmberFf14SBCanonical`):
  Applies if residue is LYS, variant is LYN, input is HZ2 or HZ3, AND
  siblings contain HZ2 and HZ3 (canonical LYN state). Predicate `true`.
  Output: input unchanged ("HZ2").

Map: 1 record. Resolver returns "HZ2". Atom's pdb_atom_name unchanged
("HZ2"). Idempotent.

Same for the atom that was HZ2 (now HZ3): canonical state, returns "HZ3"
unchanged.

### Substrate composition (later in FinalizeConstruction)

`ComposeAtomSemantic` runs. For residue 28 (variant LYN), atom HZ2:
parser produces typed identity matching LYN substrate row's HZ2 entry.
`LookupBy(LYS, 0, identity)` returns the row with PolarHKind::AmineNH,
formal_charge=0 on NZ. ✓ End-to-end works.

---

## Worked example 2: ILE δ-methyl HD1/HD2/HD3 → HD11/HD12/HD13

Fleet input has ILE residue with δ-methyl hydrogens named HD1, HD2,
HD3 (Pdb2gmx-AMBER-RTP convention). Canonical AMBER ff14SB names are
HD11, HD12, HD13.

### Pass 1 — at load

For atom HD1 in ILE residue:

```text
ctx.source              = NamingSource::Pdb2gmxAmberRtpDeviation
ctx.input_name          = "HD1"
ctx.residue_type        = AminoAcid::ILE
ctx.variant_index       = -1
ctx.terminal_state      = TerminalState::Internal
ctx.sibling_input_names = { ..., "CD", "HD1", "HD2", "HD3", ... }
                          // note: "CD" not "CD1" — also non-canonical
```

Rule iteration:

- Rule `IleCdToCD1_NonCanonicalShift` (source `AmberFf14SBCanonical`):
  Applies if residue is ILE, source is Pdb2gmxAmberRtpDeviation, input is
  "CD", AND siblings contain "CD" (not "CD1"). Predicate `false` (input
  is "HD1", not "CD"). Doesn't apply to this atom.
- Rule `IleHdMethylShift_NonCanonicalShift` (source
  `AmberFf14SBCanonical`):
  Applies if residue is ILE, source is Pdb2gmxAmberRtpDeviation, input
  is HD1/HD2/HD3, AND siblings contain HD1+HD2+HD3 but NOT HD11 (signaling
  non-canonical state — canonical state would have HD11 in siblings).
  Predicate `true`. Output: "HD11" (for input "HD1").

Map: 1 record. Resolver returns "HD11".

For atom HD2: same rule fires, output "HD12".
For atom HD3: same rule fires, output "HD13".
For atom CD: `IleCdToCD1_NonCanonicalShift` fires; output "CD1".

After Pass 1: atoms have HD11, HD12, HD13, CD1. ✓

For a HYPOTHETICAL canonical ILE input (siblings include HD11, HD12,
HD13, CD1), the same rule's predicate evaluates `false` — siblings
contain HD11, so non-canonical state predicate doesn't match. Rule
doesn't fire. Other rules don't fire either (input is canonical).
Map empty; resolver detects canonical → returns input unchanged.
Idempotent. ✓

### Pass 2 — post-protonation

ILE has no titratable variant; variant_index stays -1. Pass 2 is a
no-op for ILE atoms (no rule with predicate requiring variant_index >= 0
fires). Atom names unchanged.

---

## Worked example 3: HID variant pass-through (canonical input)

`1ubq` HIS68 protonated by Reduce as HID (Nδ1 protonated, no HE2).
Atoms include HD1; HE2 absent. cifpp parses; PdbFileReader writes
canonical names directly.

### Pass 1 — at load

For atom HD1 in HIS residue:

```text
ctx.source              = NamingSource::CifppPdbInput
ctx.input_name          = "HD1"
ctx.residue_type        = AminoAcid::HIS
ctx.variant_index       = -1
ctx.terminal_state      = TerminalState::Internal
ctx.sibling_input_names = { "N", "H", "CA", "HA", "CB", ..., "ND1", "HD1",
                            "CG", "CE1", "HE1", "CD2", "HD2", "NE2" }
                          // note: HE2 absent — HID variant
```

Rule iteration:

- Source-context rules: none with `Pdb2gmxAmberRtpDeviation` predicate fire
  (source is `CifppPdbInput`).
- Variant rules: none requiring `variant_index >= 0` fire (variant unresolved).
- Pass-through rule for canonical HIS HD1: applies if input is HD1, residue
  is HIS, siblings indicate HID-form (HD1 present, HE2 absent). Predicate
  `true`. Output: input unchanged ("HD1").

Map: 1 record. Resolver returns "HD1".

After Pass 1: HD1 unchanged. All HIS68 atoms unchanged.

### Pass 2 — post-protonation

`ResolveProtonationStates(true)` set variant_index = 0 (HID).

For HD1:

```text
ctx.variant_index = 0     // HID
ctx.input_name    = "HD1"
ctx.sibling_input_names = { ..., "ND1", "HD1", ..., "HD2", ... }
                          // post-Pass-1; same set since Pass 1 didn't rename
```

- Rule `HidCanonicalHd1PassThrough`: applies if residue is HIS, variant
  is HID, input is HD1, siblings indicate canonical HID (HD1 present, HE2
  absent). Predicate `true`. Output: input unchanged.

Map: 1 record. Returns "HD1". ✓ Idempotent.

For atom HD2 (Cδ2's hydrogen — not the variant-specific Hδ1):

- Rule `HidCanonicalHd2PassThrough`: applies if residue is HIS variant
  HID, input is HD2 (the ring CH atom). Predicate `true`. Output: "HD2".

Map: 1 record. Returns "HD2". Idempotent.

### Substrate composition

Substrate's HID variant table has rows for HD1 (PolarHKind::ImidazoleNH,
formal_charge=0 on ND1). ✓ Substrate composition succeeds.

---

## Symmetry with substrate generation

The substrate generator at `tools/topology/build_semantic_tables.cpp`
already implements this same algorithm at design time, applied to
substrate field generation:

- Multiple rule sets contribute: cifpp/CCD's `pdbx_aromatic_flag`,
  RDKit's `CIPLabeler`, RDKit's `CanonicalRankAtoms`, Markley 1998
  Table 1, Markley 1998 Figure 1, project-internal synthesis
  (`SynthesisedFor*` functions).
- Per-atom-per-field, the generator collects which sources contribute
  what value (`SemanticSourceWitness` records, up to 4 per field).
- The reconciliation precedence at
  `topology-substrate-implementation-plan-2026-05-05.md` §"Reconciliation
  precedence" is the per-case resolution method body — it's the
  documented decision policy: "RDKit CIPLabeler primary; Markley Figure
  1 cross-check; CCD pdbx_stereo_config cross-check; disagreement logs
  to RDKit-primary."
- The output (the substrate table) is the canonical truth; the witness
  stash is the audit trail (currently in `LegacyAmberSemanticTables.log.txt`).

The runtime canonicalisation architecture is the same shape, applied
at load time to atom names rather than at design time to substrate
fields. One algorithm, two scopes.

This symmetry is what makes the architecture trustworthy: it's not a
new pattern invented for naming; it's the existing substrate-generator
pattern recognised and applied at runtime.

---

## Where it lives in code

```text
src/NamingRegistry.h
src/NamingRegistry.cpp                    [refactored — replaces flat-map
                                           with NamingApplicator + rules + resolver]
src/NamingRegistry_Rules.cpp              [optional — separate file holding the
                                           rule definitions if NamingRegistry.cpp
                                           grows large]

tests/test_naming_registry.cpp            [extended — idempotency, source-aware,
                                           content-aware, conflict-resolution,
                                           fail-on-unknown tests per rule]
```

NOT touched, NOT extended:

- `src/Atom.h` — no new fields.
- `src/Protein.{h,cpp}` — no new methods. `FinalizeConstruction` continues to
  call into `NamingRegistry`/`NamingApplicator` for canonicalisation; the
  body of those calls changes shape, the signature does not.
- `src/LegacyAmberTopology.{h,cpp}` — unchanged. (The applicator may need to
  consult `AminoAcidType::atoms` for the canonicality oracle; that's a read,
  not a write to the topology.)
- All calculators — unchanged.

The applicator object is its own thing. It owns its rules, its resolver, its
diagnostics. It reads from `AminoAcidType` for the canonicality oracle. It
writes nothing onto Atom/Protein/Topology except by returning a string that
the caller assigns onto `atom.pdb_atom_name`.

---

## Open design questions for emergent resolution

These are questions the agent will encounter and may need to raise back
for design discussion. Naming them in advance so they don't surprise:

1. **Sibling-snapshot semantics.** When `ApplyResidue` snapshots sibling
   names at the start of the pass, atoms processed later in the pass see
   the original sibling set, not the partially-renamed set. This is by
   design — shift-pair rules read the non-canonical pattern from the
   original siblings to detect their applicability. Confirm: snapshot
   is per-pass-per-residue, not global, not per-atom.

2. **Canonicality oracle implementation.** The applicator needs to know
   "is this name already canonical for this residue+variant+terminal?"
   to detect idempotent inputs. Two implementation options:
   (a) Read directly from `AminoAcidType::atoms` plus the substrate's
       per-variant tables. Cleanest; substrate is the authority.
   (b) Hand-coded canonicality predicate in the applicator. Decoupled
       from substrate; more code; risk of drift.
   Recommend (a). Open: the applicator needs access to `AminoAcidType`
   and the substrate tables — pass references at construction, or read
   from globals.

3. **`std::function` per rule vs concrete `NamingRule` ABC.** Trade-off
   between concise rule definitions (function-pointer-style) and
   compile-time type safety (subclassing). Either works for the initial
   sketch; agent picks; review at codex.

4. **Diagnostic format on fail-loud.** Standard format:
   `"NamingApplicator: atom 'INPUT' in residue RES SEQ chain CHAIN
   under source SRC: <reason>. Map: <records>. Context: <context>.
   See spec/plan/naming-applicator-architecture-sketch-2026-05-06.md."`
   Tune wording during implementation.

5. **Resolution method body organisation.** As cross-rule-set conflicts
   are encountered, branches accumulate in `Resolve()`. At what point
   does the method get split into helper functions per resolution case?
   Heuristic: when the body exceeds ~50 lines, extract per-case helpers.
   Each helper unit-tested.

6. **What about rules that don't fit the per-atom shape?** Some
   canonicalisation cases might require knowing the residue's ENTIRE
   atom set as input/output (e.g., a hypothetical rule that says "if
   residue has non-canonical numbering, renumber the whole methyl group
   atomically"). For now, all observed cases are per-atom with siblings
   in context. If a residue-level rule type is needed, raise it as
   emergent issue — the architecture probably extends naturally to a
   second rule type (`NamingResidueRule` with `Applies(residue_ctx)`
   and `Output(residue_ctx) -> map<atom_index, string>`).

7. **Transient map persistence for triage.** Default: discarded after
   resolve. Debug-logging-mode: side-channel emit. Confirm default
   off; emit format TBD when first used in anger.

8. **Class naming.** The current `NamingRegistry` is a flat map; the
   new shape is more than a registry — it's registry + applicator +
   resolver. Possible renames: `NamingApplicator`, `AtomNameCanonicaliser`,
   `NamingTranslator`. Or keep `NamingRegistry` and let the architecture
   inside it carry the meaning. Default: keep `NamingRegistry` for
   minimal disruption; rename later if appropriate.

9. **`ResolveProtonationStates` interaction.** The existing protonation
   inference reads `pdb_atom_name` strings to detect variants. Per the
   user's locked decision, that path stays as-is for Phase 1. When we
   eventually revisit (post-Phase-1, separate slice), the same
   architecture applies — variant detection becomes rules with
   predicates over (residue, atom-presence-pattern), resolver picks
   variant. Out of scope for this refactor.

---

## Required tests

**Property tests on the architecture:**

- **Idempotency.** For every canonical-state input the substrate's
  `AminoAcidType::atoms` enumerates, applying the applicator returns
  the input unchanged. Iterate over all 20 standard residues + 10
  variants + cap atoms; for each canonical atom name, construct the
  context, assert `Apply(ctx) == input_name`.

- **Source-aware correctness.** A name like `HD1` in PRO under
  `Pdb2gmxAmberRtpDeviation` source returns `HD2` (shift); the same
  name under `CifppPdbInput` source returns `HD1` (canonical, no
  shift). Different source, different rule fires.

- **Content-aware correctness.** A name like `HZ2` in LYS:
  - Siblings `{HZ1, HZ2}` (non-canonical LYN) → "HZ3".
  - Siblings `{HZ2, HZ3}` (canonical LYN) → "HZ2" (idempotent).
  - Siblings `{HZ1, HZ2, HZ3}` (canonical charged LYS) → "HZ2"
    (idempotent under different rule).

- **Resolution-method correctness.** For each branch in `Resolve()`,
  construct a context that should trigger that branch and assert the
  branch's documented output. Adding a new branch = adding a test.

- **Fail-on-unknown.** Construct a context where input is non-canonical
  AND no rule applies (e.g., a residue-name + atom-name combination not
  covered by any current rule). Assert the applicator aborts via the
  project's `fprintf(stderr, "FATAL:") + std::abort()` pattern. Use
  `EXPECT_DEATH` or equivalent.

- **Conflict-resolution.** Construct a context where two rules from
  different sources both fire. Assert the resolver picks correctly
  per documented branch.

**Existing behaviour preservation:**

- All Bundle B integration tests stay green: `LegacyAmberSemanticIntegration.*`
  (8 tests).
- All existing `tests/test_naming_registry.cpp` tests stay green
  (current behaviour preserved for non-shift-pair cases).
- ctest test count rises only by the additions; no regressions.
- StringBarrier 5/5; nm RDKit count = 0.
- Generator regen byte-identical (the applicator refactor doesn't
  touch the generator).

---

## Migration plan

The refactor reshapes `NamingRegistry`'s internals while preserving
its external call sites. Steps:

1. **Inventory existing rules.** Every `AddAtomNameRule(...)` call in
   the current `NamingRegistry::InitialiseAtomNameRules` becomes a
   `NamingRule{...}` entry in the new `rules_` vector. Map each
   existing source-context (`Charmm`, `Standard`, `Amber`) to a typed
   `NamingSource` enum value that names what the rule actually
   represents (`Pdb2gmxAmberRtpDeviation`, `AmberFf14SBCanonical`,
   etc.).

2. **Rewrite predicates.** For each existing rule, define an `Applies`
   predicate that does what the old per-atom lookup did, but with
   sibling-set awareness for shift-pair rules. The 12 fleet-vetted
   "CHARMM-port" rules from Bundle B's commit `34ad91a` become rules
   with sibling-aware predicates that fire only on non-canonical
   sibling sets.

3. **Define the canonicality oracle.** A small function that returns
   the canonical chain-form atom set for a given (residue, variant,
   terminal). Reads from `AminoAcidType::atoms` plus the substrate's
   per-variant tables.

4. **Implement `Apply(ctx)` and `Resolve(map, ctx)`.** Resolve's body
   starts with the few cases observed today (single-rule firing → use
   it; canonical-input idempotency; multiple-rules conflict → branch
   by source pair with documented project decision).

5. **Replace caller sites.** `PdbFileReader.cpp:131`,
   `FullSystemReader.cpp:876` (and similar), `OrcaRunLoader.cpp:206`
   (and similar), `Protein::FinalizeConstruction`'s post-protonation
   pass currently calling `RecanonicaliseAfterProtonation` —
   all migrate to `applicator.Apply(ctx)` (or `applicator.ApplyResidue`
   for batch-residue calls).

6. **Remove `RecanonicaliseAfterProtonation`.** Its job — applying
   variant-aware rules after protonation — becomes Pass 2 of the
   uniform applicator architecture. The function deletes; its rules
   migrate to variant-aware rule entries.

7. **Add property tests** (idempotency, source-aware, content-aware,
   resolution-method, fail-on-unknown, conflict-resolution).

8. **Verify Bundle B tests stay green.** ctest 367+/367+ pass;
   StringBarrier 5/5; nm RDKit = 0; generator regen byte-identical;
   1Z9B fleet variance + ILE δ-methyl + HID/HIE/HIP variant chemistry
   all still pass.

Estimated scope: ~600-900 LoC including rule definitions, resolver
body, canonicality oracle, applicator infrastructure, and tests.
Single commit (or 2-3 atomic commits per topic). One agent dispatch.
Codex review at xhigh at the end.

---

## What's emergent vs what's locked

**Locked (do not relitigate during implementation):**

- The architecture is rule-set-preserving + per-atom transient map +
  explicit resolution method. Not first-match-wins. Not flat lookup.
  Not chained rewrites with guard checks.
- The applicator is its own thing. NOT on `Protein`, NOT on
  `LegacyAmberTopology`, NOT on `Atom`.
- Per-atom map is transient. No persistence onto Atom.
- Two passes through the applicator at `FinalizeConstruction`, with
  progressively richer context. Emergent from rule predicates, not
  architectural.
- Symmetric with substrate generator's reconciliation precedence
  pattern.
- Existing Bundle B behaviour preserved (1Z9B, 1P9J, 1ubq,
  HID/HIE/HIP variant chemistry, GLY HA2, NTERM/CTERM caps).

**Emergent (raise as user-input questions if encountered):**

- Resolution method body organisation (per-case helpers when
  body grows).
- Canonicality oracle source (substrate-direct read vs hand-coded).
- Class renaming (NamingRegistry vs NamingApplicator vs other).
- New context fields needed by future rule types.
- Residue-level rule type (if shift-pair architecture proves
  insufficient for some case).
- Debug-logging-mode emit format.
- Behavioural drift in any existing test fixture (any fixture diff
  against Bundle B baseline) — investigate before committing.

The agent's marching orders: implement the architecture as sketched,
following the data flow + worked examples + migration plan. If a
specific implementation question arises that the sketch doesn't
cover, raise it back for design discussion. Don't invent
architecture; don't paper over inconsistency; don't fall back to
chained-rewrites or first-match-wins under any circumstance.
