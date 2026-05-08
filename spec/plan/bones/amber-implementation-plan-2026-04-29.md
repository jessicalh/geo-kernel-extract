# AMBER Charge Implementation Plan

Date: 2026-04-29

Status: implementation plan for the bounded AMBER correctness slice. Six
steps, each its own commit. Steps 1–3 are pure cleanup with no NPY
drift; steps 4–6 add runtime AmberTools generation for cases the flat
ff14SB table cannot represent.

## Working-tree status (2026-04-29 evening)

```text
Step 1   GREEN in working tree   (AnalyzeFlatTableCoverage pure predicate)
Step 2   GREEN in working tree   (ResolveAmberChargeSource dispatch;
                                  all AMBER loaders routed)
Step 3   GREEN in working tree   (ChargeModelKind::AmberPreparedPrmtop;
                                  CORRECTED scope — see step 3 section)
Step 4   GREEN in working tree   (AmberPreparedChargeSource skeleton;
                                  AmberLeapInput; typed DetectDisulfides;
                                  hard re-tleap precondition asserts)
Step 5   GREEN in working tree   (tleap subprocess; PRMTOP parse;
                                  atom mapping)
Step 6   GREEN in working tree   (capping policy ACE/NME for ASH/CYM/
                                  GLH/LYN at unsupported termini)

Acceptance gate:                  62/62 tests passing across the
                                  AMBER suite + the existing pre-slice
                                  charge gate (ChargeFF14SB*, ApbsFF14SB*,
                                  PrmtopChargeTest, OrcaRunTest,
                                  FleetLoaderTest.ChargesReturned).
Full-suite ctest (368 of 373):    368 passed. 5 pre-existing failures
                                  unrelated to this slice — see below.
Master HEAD:                      f8729e3 (post-revert clean)
Working tree:                     uncommitted; per-step commit boundaries
                                  intentional, awaiting user review.
Next task (tomorrow morning):     PHASE 0 — AMBER trajectory loader,
                                  sealing AMBER-native across all five
                                  input methods.
```

### Pre-existing test failures (NOT regressions from this slice)

The full-suite ctest run on the working tree reports 5 failures. All
five are fixture/expectation mismatches unrelated to the AMBER charge
slice; they fail identically on master HEAD `f8729e3` because they do
not depend on any code this slice touches:

```text
FleetLoaderTest.HasTenFrames
FleetLoaderTest.PositionsDifferBetweenFrames
FleetLoaderTest.FullPipelineAllFrames
SmokeTest.NoDft
JobSpecE2E.FleetLibraryDirect
```

Root cause: `tests/data/fleet/1A6J_5789/poses/ensemble.json` declares
`n_poses: 1` (per master commit `b69d55c "SmokeFleet.AllPoses: trim to
1 pose; bless regenerated"`). The trim was applied to the smoke-test
fixture, but the FleetLoader tests above were not updated to match —
they still expect 10 frames. The frame-loading code is correct; the
test expectations are stale relative to the trimmed fixture.

`FleetLoaderTest.ChargesReturned` (the AMBER acceptance-gate test
named in the topology anchor) does not depend on frame count and
passes on the working tree.

A fresh session that re-runs the full ctest will see these 5
failures. **Do not attribute them to the AMBER charge slice.** They
are tracked as pre-existing fixture-expectation drift; the right
fix is a one-line fixture regeneration or a test-expectation update,
and is out of scope for this slice. PHASE 0 (AMBER trajectory loader)
will replace this CHARMM/GROMACS test fixture with an AMBER one, at
which point these tests are re-derived against the new fixture and
the drift goes away.

Audience: external AI implementors and the human reviewer.

This document supersedes any informal "tleap-it-all-the-way" framing.
The active anchor is `spec/plan/current-topology-anchor-2026-04-29.md`,
the AMBER policy is in
`spec/plan/amber-terminal-charge-generation-2026-04-29.md`, and the
contract names come from
`spec/plan/legacy-amber-implementation-brief-2026-04-29.md`.

## Goal

Every AMBER-source load of a Protein produces an authoritative
`ForceFieldChargeTable` whose rows are sourced from a documented
AmberTools artifact, or fails with a typed reason that names the
specific residue/variant that could not be supplied.

Authoritative sources, in order of preference:

```text
1. Upstream PRMTOP supplied by the load context (--orca, --mutant).
2. ff14SB flat parameter file (data/ff14sb_params.dat) when the protein's
   typed (terminal_state, variant_index, atom_name) triples are all
   covered by stock AmberTools rows.
3. PRMTOP generated at run time by tleap from typed Protein state under
   an explicit AmberPreparationPolicy (--pdb / --protonated-pdb cases
   with non-stock terminal variants).
```

There is no fourth path. There is no extractor-local chemistry. There
is no element-default fallback.

## Non-Goals

- CHARMM/GROMACS PB radii (parallel slice, parked).
- Phosphorylated, modified, selenomet, or non-protein residues.
- Calculator sweep onto `CalculatorContract<…>` (matrix work, separate).
- Splitting `ForceFieldChargeTable` into multiple typed subclasses.
- Changing `ConformationAtom::partial_charge` / `pb_radius` projection
  semantics.

## Non-Reversal Guarantees

These are properties the slice must preserve at every step.

- The `LegacyAmberTopology` + `ForceFieldChargeTable` object model from
  `spec/plan/legacy-amber-implementation-brief-2026-04-29.md`.
- The string-barrier rule: calculators do not perform residue/atom
  identity lookup by name after construction.
- The CHARMM/GROMACS quarantine: that path keeps its
  `kCompatibilityPlaceholderPbRadiusAngstrom` warning until the parallel
  slice replaces it.
- The doors for future peer topologies and naming projections
  (CHARMM topology, IUPAC projection) — nothing in this slice closes
  them.

## Architecture Summary

The "source path moves down into the logic" discipline:

```text
Loaders populate typed state:
    Protein (residues, terminal_state, protonation_variant_index, bonds)
    ProteinBuildContext (force_field intent, upstream prmtop_path,
        protonation_tool)

A single resolver function decides which ChargeSource to construct:
    ResolveAmberChargeSource(protein, build_context, config)
    -> unique_ptr<ChargeSource>

The decision is a typed predicate over Protein state, not a loader switch.

The output ChargeSource carries its origin via ChargeModelKind, and the
resulting ForceFieldChargeTable preserves that kind plus a typed
description string for provenance.
```

New entities (additions to existing object model):

```text
src/AmberChargeResolver.{h,cpp}
    AmberSourceConfig                  config struct (paths, policy)
    AmberFlatTableCoverageVerdict      typed sum (Satisfiable | UnsupportedTerminalVariant | UnsupportedResidue | MissingAtomName)
    AnalyzeFlatTableCoverage(...)      pure predicate
    ResolveAmberChargeSource(...)      single dispatch function

src/AmberPreparedChargeSource.{h,cpp}
    AmberPreparationPolicy             enum
    AmberPreparedChargeSource          ChargeSource subclass: tleap-driven

src/AmberLeapInput.{h,cpp}             (private to the prepared source)
    GenerateAmberPdb(protein, conf, policy, work_dir)
    GenerateLeapScript(protein, policy, paths)
    ResidueMap                         extractor↔prmtop residue index map
```

Extended (additions to existing entities, no removal):

```text
ChargeModelKind
    + AmberPreparedPrmtop
    + ChargeModelKindName helper

ProteinBuildContext
    (NO CHARGE FIELDS ADDED. See Step 3 "CORRECTED 2026-04-29"
     section below — ForceFieldChargeTable is the typed authority
     on charge-source identity via Kind() and SourceDescription();
     duplicating those facts on BuildContext was rejected as a
     second-authority anti-pattern.)
```

`ChargeAssignmentResult` does not change. It remains a thin projection
over `ForceFieldChargeTable`.

## Step Plan

Each step is one commit. Steps 1–3 produce identical NPY output (no
drift). Steps 4–6 introduce numerical drift only for proteins that
previously hit the loud-fail path (i.e., proteins that did not produce
output at all). NPY parity for proteins on the existing flat-table path
is preserved through step 6.

---

### Step 1: Extract `AnalyzeFlatTableCoverage` as a pure predicate

**Goal.** Pull the verdict logic out of `ParamFileChargeSource::LoadCharges`
so the resolver can decide *before* calling LoadCharges. No behavior
change.

**Files created:**

```text
src/AmberChargeResolver.h       (interface)
src/AmberChargeResolver.cpp     (implementation; resolver added in step 2)
```

**Files modified:**

```text
src/ChargeSource.cpp
    ParamFileChargeSource::LoadCharges keeps the same external behavior
    but delegates verdict computation to AnalyzeFlatTableCoverage.

CMakeLists.txt
    add AmberChargeResolver.{h,cpp} to the library target.
```

**API shape:**

```cpp
namespace nmr {

enum class AmberFlatTableCoverageKind {
    Satisfiable,
    UnsupportedTerminalVariant,   // (terminal, ff_resname) prefix has no rows
    UnsupportedResidue,           // ff_resname has no rows in any terminal state
    MissingAtomName               // (terminal, ff_resname) has rows but atom_name missing
};

struct AmberFlatTableCoverageVerdict {
    AmberFlatTableCoverageKind kind = AmberFlatTableCoverageKind::Satisfiable;
    // Diagnostic fields, populated when kind != Satisfiable:
    std::string terminal_token;
    std::string ff_residue_name;
    std::string atom_name;
    int residue_sequence_number = 0;
    std::string chain_id;
    std::string detail;           // human-readable summary
    bool Ok() const { return kind == AmberFlatTableCoverageKind::Satisfiable; }
};

AmberFlatTableCoverageVerdict AnalyzeFlatTableCoverage(
    const Protein& protein,
    const std::string& flat_table_path);

}  // namespace nmr
```

**Behavior.**

The predicate parses the flat table once into the same row-key set
`ParamFileChargeSource::LoadFf14sbParamFile` builds. It walks every
residue and every atom of the protein, computes
`(TerminalStateToken(res.terminal_state), Ff14sbVariantResidueName(...),
 atom_candidate_in {pdb_atom_name, "H1" if NTERM and pdb_atom_name in
 {"H","HN"}})` and confirms a row exists. The first missing case sets
the verdict and returns. A satisfiable verdict means
`ParamFileChargeSource::LoadCharges` will succeed end-to-end.

`ParamFileChargeSource::LoadCharges` is refactored to call this first.
On non-satisfiable verdicts it sets the same `error_out` it sets today
(the resolver in step 2 inspects the verdict directly; the refactored
LoadCharges error path is preserved for callers that still take the
direct-source path).

**Test plan.**

```text
tests/test_amber_charge_resolver.cpp (new)
    Satisfiable on 1UBQ protonated PDB.
    UnsupportedTerminalVariant on a synthetic protein with NTERM ASH.
    UnsupportedTerminalVariant on a synthetic protein with NTERM TYM.
    UnsupportedResidue on a protein with a non-standard residue
        (synthesised; out-of-band for production but verifies the kind).
    MissingAtomName when a custom dat file has the (terminal, resname)
        but lacks one of the residue's atoms.
    Verdict diagnostic fields populated correctly.
```

**Acceptance.**

```text
cmake --build build -j2
ctest -R AmberChargeResolver
ctest -R ChargeFF14SB
ctest -R ApbsFF14SB
ctest -R PrmtopCharge
ctest -R OrcaRun
ctest -R FleetLoader.ChargesReturned

NPY parity: identical to current HEAD for blessed proteins.
```

---

### Step 2: Land `ResolveAmberChargeSource` and route AMBER loaders through it

**Goal.** The single dispatch function. Loaders stop constructing
`ChargeSource` subclasses directly; they call the resolver. Same
behavior, cleaner authority.

**Files created:** none (continues `AmberChargeResolver.cpp`).

**Files modified:**

```text
src/PdbFileReader.cpp
    BuildFromPdb and BuildFromProtonatedPdb construct a ProteinBuildContext
    populated with force_field/protonation_tool/etc., then call the
    resolver instead of constructing ParamFileChargeSource directly.

src/OrcaRunLoader.cpp
    BuildFromOrca constructs the same context, then calls the resolver.
    The upstream-PRMTOP branch in the resolver short-circuits to
    PrmtopChargeSource without further work.

src/AmberChargeResolver.{h,cpp}
    Adds AmberSourceConfig and ResolveAmberChargeSource.
```

**API shape:**

```cpp
namespace nmr {

enum class AmberPreparationPolicy {
    UseStockTermini,                                     // step 4 stub
    UseCappedFragmentsForUnsupportedTerminalVariants,    // step 6
    FailOnUnsupportedTerminalVariants                    // default until step 4 lands runtime tleap
};

struct AmberSourceConfig {
    std::string flat_table_path;                          // e.g. data/ff14sb_params.dat
    AmberPreparationPolicy preparation_policy =
        AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
    std::string tleap_path;                               // resolved from RuntimeEnvironment in step 4+
    std::string work_dir;                                 // tmp root for AmberPreparedChargeSource (step 4+)
};

std::unique_ptr<ChargeSource> ResolveAmberChargeSource(
    const Protein& protein,
    const ProteinBuildContext& build_context,
    const AmberSourceConfig& config,
    std::string& error_out);

}  // namespace nmr
```

**Resolver logic (this is the only place the dispatch lives):**

```cpp
std::unique_ptr<ChargeSource> ResolveAmberChargeSource(
        const Protein& protein,
        const ProteinBuildContext& build_context,
        const AmberSourceConfig& config,
        std::string& error_out) {

    // Branch 1: upstream PRMTOP wins, always. Never re-tleap.
    if (!build_context.prmtop_path.empty()) {
        if (!std::filesystem::exists(build_context.prmtop_path)) {
            error_out = "build_context.prmtop_path does not exist: " +
                        build_context.prmtop_path;
            return nullptr;
        }
        return std::make_unique<PrmtopChargeSource>(
            build_context.prmtop_path, ForceField::Amber_ff14SB);
    }

    // Branch 2: typed predicate over Protein state.
    auto verdict = AnalyzeFlatTableCoverage(protein, config.flat_table_path);
    if (verdict.Ok()) {
        return std::make_unique<ParamFileChargeSource>(config.flat_table_path);
    }

    // Branch 3: prepared topology via tleap.
    if (config.preparation_policy ==
            AmberPreparationPolicy::FailOnUnsupportedTerminalVariants) {
        error_out = "ff14SB flat table cannot represent this protein and "
                    "policy is FailOnUnsupportedTerminalVariants: " +
                    verdict.detail;
        return nullptr;
    }
    return std::make_unique<AmberPreparedChargeSource>(
        protein, config.preparation_policy, verdict, config);
}
```

(`AmberPreparedChargeSource` is a forward declaration in step 2; step 4
makes it real. In step 2 the ABC line exists but instantiation paths
fail loudly because policy default is FailOnUnsupportedTerminalVariants.)

**Loader call shape:**

```cpp
// Before (in PdbFileReader.cpp):
result.charges = std::make_unique<ParamFileChargeSource>(
    RuntimeEnvironment::Ff14sbParamsPath());

// After:
AmberSourceConfig source_config;
source_config.flat_table_path = RuntimeEnvironment::Ff14sbParamsPath();
source_config.preparation_policy =
    AmberPreparationPolicy::FailOnUnsupportedTerminalVariants;
std::string charge_err;
result.charges = ResolveAmberChargeSource(
    *result.protein, *result.protein->BuildContextPtr(),
    source_config, charge_err);
if (!result.charges) {
    result.error = charge_err;
    return result;
}
```

**Test plan.**

```text
tests/test_amber_charge_resolver.cpp (extended)
    ResolverPicksPrmtopWhenBuildContextHasIt
        Synthetic Protein + BuildContext with prmtop_path → PrmtopChargeSource
    ResolverPicksFlatTableForStockProtein
        1UBQ protonated → ParamFileChargeSource
    ResolverFailsForUnsupportedTerminalVariantUnderFailPolicy
        Synthetic NTERM ASH → returns nullptr with named error.
    ResolverProducesPreparedSourceUnderUseCappedFragmentsPolicy
        (Step 4 makes this test produce a non-null source; in step 2 it
         fails until policy default is moved.)

Existing tests:
    BuildFromPdb / BuildFromProtonatedPdb / BuildFromOrca paths produce
    identical ForceFieldChargeTable contents to step 1 (NPY parity).
```

**Acceptance.**

```text
cmake --build build -j2
ctest -R AmberChargeResolver
ctest -R ChargeFF14SB
ctest -R ApbsFF14SB
ctest -R PrmtopCharge
ctest -R OrcaRun
ctest -R FleetLoader.ChargesReturned

NPY parity: identical to current HEAD for blessed proteins.
```

---

### Step 3: `ChargeModelKind` extension — CORRECTED 2026-04-29

**Goal.** Add the `AmberPreparedPrmtop` discriminator to the typed charge
authority that already exists. Charge-source provenance lives on
`ForceFieldChargeTable` (already has `Kind()` and `SourceDescription()`);
this step only extends the enum so step 4–6 sources can identify
themselves.

**Scope correction (2026-04-29).** An earlier draft of this step
proposed adding `charge_source_kind` / `charge_source_description` /
`preparation_policy_name` / `preparation_reason` fields to
`ProteinBuildContext`. That was wrong:

```text
ProteinBuildContext is the protein-level provenance record: PDB source,
organism, protonation_tool, force_field NAME, prmtop_path,
tleap_script_path. Charge-source identity is not a protein-level fact;
it is the prepared charge object's typed identity. ForceFieldChargeTable
is already the authority — it carries source_force_field_, kind_,
source_description_ on construction. Adding parallel fields to
ProteinBuildContext would make it a second authority on the same
information, exactly the anti-pattern the LegacyAmberTopology landing
fixed for topology.
```

The chain that survives: `loader → resolver → ForceFieldChargeTable`.
`Protein::ForceFieldCharges().Kind()` and `.SourceDescription()` are
the canonical reads for methods text and H5 attribute emission.
`ProteinBuildContext` does not change in this step.

**Files modified:**

```text
src/ChargeSource.h
    enum class ChargeModelKind: add AmberPreparedPrmtop.
    ChargeModelKindName helper.
```

**Files NOT modified by this step:**

```text
src/ProteinBuildContext.h     (no charge fields added)
src/Protein.cpp                (PrepareForceFieldCharges does not write
                                 to build_context_; ForceFieldChargeTable
                                 already records the source identity at
                                 construction time)
```

**API shape:**

```cpp
enum class ChargeModelKind {
    AmberPrmtop,            // upstream PRMTOP via PrmtopChargeSource
    AmberPreparedPrmtop,    // PRMTOP we generated via tleap (step 4+)
    GromacsTpr,             // CHARMM/GROMACS, currently placeholder PB radii
    Ff14SBParamFile,        // flat dat file
    Preloaded,
    Unknown
};

inline const char* ChargeModelKindName(ChargeModelKind k);
```

`AmberPreparedChargeSource` (step 4) constructs itself with
`Kind() == AmberPreparedPrmtop`. Its `Describe()` includes the policy
name and verdict reason in human-readable form, which
`ForceFieldChargeTable::source_description_` carries forward. If a
later slice needs *structured* policy/reason access after the source
goes out of scope, the right place to land that is on
`ForceFieldChargeTable` (extending the prepared object), not on
`ProteinBuildContext`.

**Test plan.**

```text
tests/test_amber_charge_resolver.cpp
    AmberPreparationPolicyName.AllValuesNamed (already added in step 1).
    No new BuildContext-charge tests; that surface does not exist by
    design.
```

**Acceptance.**

```text
cmake --build build -j2
ctest -R AmberChargeResolver
ctest -R ChargeFF14SB
ctest -R ApbsFF14SB
ctest -R PrmtopCharge
ctest -R OrcaRun
ctest -R FleetLoader.ChargesReturned

NPY parity: identical to current HEAD for blessed proteins.
```

---

### Step 4: `AmberPreparedChargeSource` skeleton — script + PDB generation, no tleap call

> **Post-step-5 evolution note (2026-04-29):** Step 4 *as recorded
> below* lands the deterministic scaffolding only. **Steps 5 and 6
> have since landed on top of it**: step 5 wired the actual `tleap`
> subprocess + PRMTOP parse + atom mapping, step 6 added the capping
> policy. `AmberPreparedChargeSource::LoadCharges` is now the full
> pipeline, not a placeholder. The step-4 description below is the
> historical commit-boundary record; current behaviour is described
> in steps 5 and 6.

**Goal.** All the deterministic, testable parts of the runtime AMBER
preparation, with the actual `tleap` invocation stubbed. Tests verify
that for a given Protein + policy, the LEaP script and the generated
PDB are byte-deterministic and contain what they should.

**Files created:**

```text
src/AmberPreparedChargeSource.h
src/AmberPreparedChargeSource.cpp
src/AmberLeapInput.h
src/AmberLeapInput.cpp
```

**Files modified:**

```text
CMakeLists.txt
    add the new files.
src/AmberChargeResolver.cpp
    can now construct AmberPreparedChargeSource for non-Fail policies.
```

**Pipeline split:**

```cpp
// AmberLeapInput.h
namespace nmr::amber_leap {

struct ResidueAmberMapping {
    // Per generated-PRMTOP residue:
    //   extractor_index = NONE_FOR_CAP for inserted ACE/NME/NHE residues
    //   extractor_index = N for residues that correspond to
    //                       protein.ResidueAt(N)
    static constexpr size_t NONE_FOR_CAP = std::numeric_limits<size_t>::max();
    std::vector<size_t> extractor_index_for_prmtop_residue;  // 0-indexed parallel
};

struct LeapScriptInputs {
    std::string pdb_path;
    std::string prmtop_path;
    std::string inpcrd_path;
    std::vector<std::pair<size_t,size_t>> disulfide_residue_pairs_1based;
};

// Pure functions; no I/O of their own (they take std::ostream& or
// std::string& out so tests can capture deterministic output).

void GenerateAmberPdb(const Protein& protein,
                      const ProteinConformation& conf,
                      AmberPreparationPolicy policy,
                      const AmberFlatTableCoverageVerdict& verdict,
                      std::ostream& pdb_out,
                      ResidueAmberMapping& map_out);

void GenerateLeapScript(const LeapScriptInputs& inputs,
                        std::ostream& script_out);

}  // namespace nmr::amber_leap
```

**Generated PDB rules.**

```text
- Atom records emitted in extractor atom-index order, except where caps
  are inserted.
- Each residue's residue_name is the AMBER unit name derived from
  (residue.type, residue.protonation_variant_index, residue.terminal_state)
  using the typed AminoAcidType variants table. HIS becomes HID/HIE/HIP,
  CYS becomes CYS or CYX (CYX iff a Disulfide bond from CovalentTopology
  touches its SG), ASP becomes ASP or ASH, etc.
- ICODE preserved.
- TER records between chains, derived from Residue::chain_id.
- Cap residues (ACE/NME/NHE) inserted only when policy allows AND verdict
  named the chain end as unsupported. Cap residues get atom names from
  the AMBER libraries (ACE: HH31 HH32 HH33 CH3 C O; NME: N H CH3 HH31
  HH32 HH33; NHE: N H1 H2). Cap atoms have synthesized positions from
  internal AMBER ideal geometry (this slice writes ideal-geometry caps;
  refinement is upstream's job).
- atom_name on regular atoms is the residue's existing pdb_atom_name
  (already AMBER convention on AMBER-source loads — this slice does not
  cross the CHARMM line).
- ResidueAmberMapping records, for each PDB residue index in emission
  order, which extractor residue it corresponds to (or NONE_FOR_CAP).
```

**Generated LEaP script rules.**

```text
source leaprc.protein.ff14SB
set default PBRadii mbondi2
mol = loadPdb <pdb_path>
<for each disulfide pair (ri, rj) in 1-based PRMTOP indexing>
bond mol.<ri>.SG mol.<rj>.SG
saveamberparm mol <prmtop_path> <inpcrd_path>
quit
```

The disulfide pair list comes from `protein.LegacyAmber().Bonds()`
filtered by `BondCategory::Disulfide`, with 1-based PRMTOP indices
derived from the ResidueAmberMapping (cap residues shift the indexing).

**`AmberPreparedChargeSource` skeleton:**

```cpp
class AmberPreparedChargeSource : public ChargeSource {
public:
    AmberPreparedChargeSource(const Protein& protein,
                              AmberPreparationPolicy policy,
                              AmberFlatTableCoverageVerdict reason,
                              const AmberSourceConfig& config);

    ForceField SourceForceField() const override { return ForceField::Amber_ff14SB; }
    ChargeModelKind Kind() const override { return ChargeModelKind::AmberPreparedPrmtop; }
    std::string Describe() const override;

    std::vector<AtomChargeRadius> LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const override;

    // For step 4, also exposed for testing without tleap:
    std::string GeneratedPdb(const ProteinConformation& conf) const;
    std::string GeneratedLeapScript(const std::string& pdb_path,
                                    const std::string& prmtop_path,
                                    const std::string& inpcrd_path) const;
    const amber_leap::ResidueAmberMapping& ResidueMapping() const;

private:
    // ... see step 5 for the LoadCharges body.
};
```

In step 4, `LoadCharges` checks the hard preconditions and returns a
named "tleap invocation not yet implemented" error. The tests exercise
`GeneratedPdb` and `GeneratedLeapScript` directly.

**Hard preconditions enforced in LoadCharges:**

```cpp
if (protein.HasForceFieldCharges()) {
    fprintf(stderr, "FATAL: AmberPreparedChargeSource on Protein that "
                    "already has charges (re-tleap protection).\n");
    std::abort();
}
if (!protein.BuildContext().prmtop_path.empty()) {
    fprintf(stderr, "FATAL: AmberPreparedChargeSource when build context "
                    "names an upstream PRMTOP (re-tleap protection).\n");
    std::abort();
}
if (config_.tleap_path.empty()) {
    error_out = "AmberPreparedChargeSource: tleap_path is empty (set "
                "RuntimeEnvironment AMBERHOME or pass --tleap)";
    return {};
}
```

**Test plan.**

```text
tests/test_amber_leap_input.cpp (new)
    GeneratedPdb_StockProtein_ByteEqualToFixture
        1UBQ protonated → known PDB string in tests/data/amber/.
    GeneratedPdb_HisVariantUsesHidHieHipNames
    GeneratedPdb_DisulfidesEmitCYX
    GeneratedPdb_TerRecordsBetweenChains
    GeneratedPdb_NTermAshUnderCappingPolicyInsertsAce
    GeneratedPdb_CTermAshUnderCappingPolicyInsertsNme
    GeneratedPdb_FailPolicyDoesNotInsertCaps                // emits as-is, tleap will fail; this slice doesn't run tleap
    ResidueAmberMapping_NoCaps_OneToOne
    ResidueAmberMapping_WithCap_ExtractorIndexShifted

    GeneratedLeapScript_NoDisulfides_ProducesExpectedString
    GeneratedLeapScript_WithDisulfides_EmitsBondCommands
    GeneratedLeapScript_DisulfideIndicesAccountForCaps      // ri, rj are 1-based PRMTOP indices

tests/test_amber_charge_resolver.cpp (extended)
    HardPrecondition_ReTleapAbortsWhenChargesAlreadyPresent
        Death test: protein with ForceFieldChargeTable already attached.
    HardPrecondition_ReTleapAbortsWhenBuildContextHasPrmtop
        Death test: build_context.prmtop_path non-empty.
    LoadChargesReturnsNamedErrorWhenTleapPathEmpty
```

**Acceptance.**

```text
cmake --build build -j2
ctest -R AmberLeapInput
ctest -R AmberChargeResolver
ctest -R ChargeFF14SB
ctest -R ApbsFF14SB
ctest -R PrmtopCharge
ctest -R OrcaRun
ctest -R FleetLoader.ChargesReturned

NPY parity: identical to current HEAD for blessed proteins.
(No protein in the blessed set hits the prepared-source path because
default policy is FailOnUnsupportedTerminalVariants.)
```

---

### Step 5: Tleap invocation, PRMTOP parse, atom mapping

**Goal.** End-to-end runtime AMBER preparation. Wired but only fires for
proteins whose verdict is non-satisfiable AND policy is non-Fail.

**Files modified:**

```text
src/AmberPreparedChargeSource.cpp
    LoadCharges body landed.
src/RuntimeEnvironment.{h,cpp}
    + std::string TleapPath() const;       // resolved at startup (e.g. via AMBERHOME)
    + std::string AmberPreparedWorkDir() const;
```

**LoadCharges body:**

```cpp
std::vector<AtomChargeRadius> AmberPreparedChargeSource::LoadCharges(
        const Protein& protein,
        const ProteinConformation& conf,
        std::string& error_out) const {

    OperationLog::Scope scope("AmberPreparedChargeSource::LoadCharges",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " policy=" + AmberPreparationPolicyName(policy_) +
        " reason=" + reason_.detail);

    // 1. Hard preconditions (death asserts as listed above).
    // 2. Make a deterministic work directory inside config_.work_dir.
    std::string work_dir = MakeWorkDir(protein, conf);
    fs::create_directories(work_dir);

    std::string pdb_path     = work_dir + "/input.pdb";
    std::string script_path  = work_dir + "/leap.in";
    std::string prmtop_path  = work_dir + "/prep.prmtop";
    std::string inpcrd_path  = work_dir + "/prep.inpcrd";
    std::string log_path     = work_dir + "/tleap.log";

    // 3. Generate PDB and LEaP script. Save mapping for atom mapping later.
    {
        std::ofstream pdb_out(pdb_path);
        if (!pdb_out) { error_out = "cannot write " + pdb_path; return {}; }
        amber_leap::GenerateAmberPdb(protein, conf, policy_, reason_,
                                     pdb_out, residue_mapping_);
    }
    {
        std::ofstream script_out(script_path);
        amber_leap::LeapScriptInputs inputs;
        inputs.pdb_path = pdb_path;
        inputs.prmtop_path = prmtop_path;
        inputs.inpcrd_path = inpcrd_path;
        inputs.disulfide_residue_pairs_1based =
            DisulfidePairsFromTopology(protein, residue_mapping_);
        amber_leap::GenerateLeapScript(inputs, script_out);
    }

    // 4. Run tleap.
    std::string cmd = config_.tleap_path + " -f " + script_path +
                      " > " + log_path + " 2>&1";
    int rc = std::system(cmd.c_str());
    if (rc != 0 || !fs::exists(prmtop_path)) {
        error_out = "tleap failed (rc=" + std::to_string(rc) +
                    "); see " + log_path;
        return {};
    }

    // 5. Parse PRMTOP CHARGE / RADII / ATOM_NAME / RESIDUE_LABEL /
    //    RESIDUE_POINTER. (Reuse the parser from PrmtopChargeSource;
    //    extract its ReadPrmtopSection / new ReadPrmtopStrings to a
    //    shared helper if needed.)
    auto prmtop = ReadPrmtopFields(prmtop_path);
    if (!prmtop.Ok()) {
        error_out = "PRMTOP parse failed: " + prmtop.error;
        return {};
    }

    // 6. Atom mapping.
    return MapPrmtopToExtractor(protein, prmtop, residue_mapping_, error_out);
}
```

**Atom-mapping algorithm:**

```cpp
std::vector<AtomChargeRadius> MapPrmtopToExtractor(
        const Protein& protein,
        const PrmtopFields& prmtop,
        const amber_leap::ResidueAmberMapping& mapping,
        std::string& error_out) {

    const size_t n_extractor_atoms = protein.AtomCount();
    std::vector<AtomChargeRadius> result(n_extractor_atoms);
    std::vector<bool> seen(n_extractor_atoms, false);

    // For each extractor residue, build an index of its (atom_name -> atom_index).
    std::vector<std::map<std::string, size_t>> name_index_by_extractor_residue(
        protein.ResidueCount());
    for (size_t ai = 0; ai < n_extractor_atoms; ++ai) {
        const Atom& atom = protein.AtomAt(ai);
        const Residue& res = protein.ResidueAt(atom.residue_index);
        auto& m = name_index_by_extractor_residue[atom.residue_index];
        if (m.count(atom.pdb_atom_name) > 0) {
            error_out = "extractor protein has duplicate atom name " +
                        atom.pdb_atom_name + " in residue " +
                        std::to_string(res.sequence_number);
            return {};
        }
        m[atom.pdb_atom_name] = ai;
    }

    // PRMTOP residue 0 = mapping.extractor_index_for_prmtop_residue[0].
    // Walk PRMTOP atoms.
    for (size_t prmtop_ai = 0; prmtop_ai < prmtop.atom_names.size(); ++prmtop_ai) {
        size_t prmtop_ri = prmtop.residue_index_for_atom[prmtop_ai];
        size_t ext_ri = mapping.extractor_index_for_prmtop_residue[prmtop_ri];
        if (ext_ri == amber_leap::ResidueAmberMapping::NONE_FOR_CAP) {
            continue;  // ACE/NME/NHE inserted-cap atom, ignore.
        }

        const std::string& aname = prmtop.atom_names[prmtop_ai];
        auto& m = name_index_by_extractor_residue[ext_ri];
        auto it = m.find(aname);
        if (it == m.end()) {
            error_out = "PRMTOP atom " + aname + " in residue " +
                        prmtop.residue_labels[prmtop_ri] + " (extractor index " +
                        std::to_string(ext_ri) +
                        ") has no match in extractor protein";
            return {};
        }
        size_t ext_ai = it->second;
        if (seen[ext_ai]) {
            error_out = "PRMTOP atom " + aname + " mapped twice to extractor "
                        "atom index " + std::to_string(ext_ai);
            return {};
        }
        seen[ext_ai] = true;
        result[ext_ai].partial_charge =
            prmtop.charges[prmtop_ai] / 18.2223;
        result[ext_ai].pb_radius = prmtop.radii[prmtop_ai];
        result[ext_ai].status = ChargeAssignmentStatus::Matched;
    }

    // Every extractor atom must have been seen.
    for (size_t ai = 0; ai < n_extractor_atoms; ++ai) {
        if (!seen[ai]) {
            const Atom& atom = protein.AtomAt(ai);
            const Residue& res = protein.ResidueAt(atom.residue_index);
            error_out = "extractor atom " + atom.pdb_atom_name +
                        " in residue " +
                        std::to_string(res.sequence_number) +
                        " has no PRMTOP counterpart";
            return {};
        }
    }

    return result;
}
```

**Provenance write-back.** After successful LoadCharges,
`Protein::PrepareForceFieldCharges` populates the BuildContext fields
including `preparation_policy_name = AmberPreparationPolicyName(policy)`
and `preparation_reason = reason.detail`.

**Test plan.**

```text
tests/test_amber_prepared_charge_source.cpp (new; tagged
                                              AMBER_REQUIRES_TLEAP)
    These tests are gated by RuntimeEnvironment::TleapPath() being
    available; otherwise they SKIP.

    Roundtrip_StockProtein_AgreesWithFlatTable
        Take 1UBQ, force AmberPreparedChargeSource (bypassing the
        resolver's flat-table preference), verify resulting charges
        agree with the flat-table charges to 1e-6 per atom.

    Roundtrip_NTermHisVariants
        Synthetic NTERM HID/HIE/HIP proteins → expected ff14SB
        terminal HID/HIE/HIP charges.

    Roundtrip_Disulfide
        Synthetic two-CYS protein with a CovalentTopology Disulfide bond
        → both SGs come back as CYX charges, not CYS.

    AtomMapping_FailsOnUnmatchedExtractorAtom
        Synthetic protein with a stray atom (constructed mid-test) →
        loud error naming the unmatched atom.

    AtomMapping_DropsCapAtomsWithoutLeak
        N-terminal ASH under capping policy → result has charges only
        for the original residue's atoms; ACE atoms not present in the
        result vector and don't shift indices.
```

**Acceptance.**

```text
cmake --build build -j2
ctest -R AmberPreparedChargeSource
ctest -R AmberLeapInput
ctest -R AmberChargeResolver
ctest -R ChargeFF14SB
ctest -R ApbsFF14SB
ctest -R PrmtopCharge
ctest -R OrcaRun
ctest -R FleetLoader.ChargesReturned

NPY parity: identical to current HEAD for blessed proteins (the existing
blessed set has no NTERM ASH/CYM/GLH/LYN, so the prepared-source path is
not hit).
```

---

### Step 6: Capping policy as a real, tested policy

**Goal.** `UseCappedFragmentsForUnsupportedTerminalVariants` produces
correct charges for proteins with unsupported terminal protonation
variants, with the methods-level decision recorded as provenance.

**Files modified:**

```text
src/AmberLeapInput.cpp
    GenerateAmberPdb under UseCappedFragmentsForUnsupportedTerminalVariants:
        For each chain end whose verdict named that residue as
        UnsupportedTerminalVariant AND whose variant is one of
        {ASH, CYM, GLH, LYN}:
            Emit ACE before the residue (for NTERM cases) or NME after
            (for CTERM cases). Cap residue atoms get ideal AMBER geometry
            positions. ACE/NME residue indices increment the PRMTOP-side
            counter; ResidueAmberMapping marks them NONE_FOR_CAP.
            The original residue is emitted with its INTERNAL variant
            name (ASH/CYM/GLH/LYN), not the unsupported terminal one.

src/AmberPreparedChargeSource.cpp
    Describe() includes the policy name and a list of the capped
    chain ends.

src/AmberChargeResolver.cpp
    No change; the resolver already routes based on policy.
```

**Cap atom geometry.** For step 6, ideal-geometry caps are sufficient:

```text
ACE cap at N-terminus of residue R:
    placed along the -N-CA backbone direction, opposite the rest of
    the chain
    ACE.C bonded to R.N, ACE.O double-bonded to ACE.C
    ACE.CH3 with three Hs in tetrahedral geometry
    coordinates: linearly extrapolated from R.N - R.CA direction with
    standard AMBER bond lengths (CH3-C 1.522, C-O 1.229, C-N 1.335)
    and tetrahedral angles for CH3 hydrogens.

NME cap at C-terminus of residue R:
    placed along the C-N direction, ahead of the chain
    NME.N bonded to R.C, NME.H on N, NME.CH3 with three Hs.

NHE cap at C-terminus of residue R:
    NHE.N bonded to R.C with two Hs (no methyl).
```

These are AMBER-standard geometries. tleap accepts them and recomputes H
positions if needed. We do not refine — refinement is upstream's job.

Default chosen for `UseCappedFragmentsForUnsupportedTerminalVariants`:
ACE for NTERM, NME for CTERM. NHE is not used by default (NME is more
common in published methods text).

**Provenance.**

```text
build_context.preparation_policy_name = "UseCappedFragmentsForUnsupportedTerminalVariants"
build_context.preparation_reason = "NTERM ASH at residue 1 (chain A) capped with ACE; CTERM GLH at residue 76 (chain A) capped with NME"
```

The reason string lists every cap inserted. Methods text reads from it
verbatim.

**Test plan.**

```text
tests/test_amber_leap_input.cpp (extended)
    GeneratedPdb_NTermAshWithCappingHasAceWithIdealGeometry
        ACE atoms placed correctly; bond lengths/angles within tolerance.
    GeneratedPdb_CTermAshWithCappingHasNme
    GeneratedPdb_NTermAndCTermBothCappedSimultaneously
    GeneratedPdb_CappingDoesNotApplyToHidHieHipNterm
        These are stock-supported; capping policy must NOT cap them.

tests/test_amber_prepared_charge_source.cpp (extended; AMBER_REQUIRES_TLEAP)
    Roundtrip_NTermAshWithCapping_ProducesCorrectInternalAshCharges
        Charges on the original residue's atoms match
        flat-table INTERNAL ASH within 1e-6.
    Roundtrip_CTermGlhWithCapping_ProducesCorrectInternalGlhCharges
    BuildContextRecordsCappingDecision
    AtomMapping_CapAtomsExcludedFromResultVector
```

**Acceptance.**

```text
cmake --build build -j2
ctest -R AmberPreparedChargeSource
ctest -R AmberLeapInput
ctest -R AmberChargeResolver
ctest -R ChargeFF14SB
ctest -R ApbsFF14SB
ctest -R PrmtopCharge
ctest -R OrcaRun
ctest -R FleetLoader.ChargesReturned

NPY parity: identical to current HEAD for blessed proteins.

Curated audit (not gated; run by user before declaring this slice done):
    For every protein in tests/data/curated/ (or the calibration set if
    accessible), run BuildFromProtonatedPdb under both policies:
        FailOnUnsupportedTerminalVariants  (current default)
        UseCappedFragmentsForUnsupportedTerminalVariants
    Record:
        - count of proteins requiring capping
        - per-protein list of chain ends capped
        - per-protein delta in net_charge (should be 0 — caps are neutral)
        - per-protein NPY hash for stock proteins (should be identical
          to flat-table baseline)
```

---

## Cross-Cutting Concerns

### Re-tleap protection

Two layers:

1. **Resolver dispatch**: branch 1 (upstream PRMTOP) returns
   `PrmtopChargeSource` and never reaches branch 3. The resolver is the
   only place AMBER ChargeSources are constructed.
2. **Source preconditions**: `AmberPreparedChargeSource::LoadCharges`
   asserts that the Protein has no charges and BuildContext has no
   prmtop_path. Violation → `std::abort()`.

The two layers are redundant by intent. The resolver is the
contract; the asserts are a structural guard against a future caller
that bypasses the resolver.

### Atom mapping correctness

The single point of truth is `ResidueAmberMapping`, built during PDB
generation. After tleap, the PRMTOP's residue order matches the
generated PDB's residue order (tleap's `loadPdb` preserves order). The
mapping's `extractor_index_for_prmtop_residue` array gives every
PRMTOP residue its extractor counterpart or NONE_FOR_CAP.

Failure modes that abort:

```text
extractor atom unmatched               → loud error naming the atom
non-cap PRMTOP atom unmatched          → loud error naming the atom
two-PRMTOP-atoms-to-one-extractor-atom → loud error
extractor protein has duplicate atom name in same residue → loud error (this is a Protein-construction bug, not an AMBER issue)
```

### Multi-chain proteins

The generated PDB emits TER records between chains. tleap's loadPdb
respects them and creates separate "molecules" within the unit; saveamberparm
combines them in PRMTOP atom order matching PDB order. ResidueAmberMapping
is constructed in the same order, so the mapping is consistent.

### Insertion codes (ICODE)

Preserved in the generated PDB. PRMTOP RESIDUE_LABEL doesn't carry ICODE,
but residue order does, so the index-based mapping is unaffected.

### Disulfide detection

Disulfide pairs come from `protein.LegacyAmber().Bonds()` filtered by
`BondCategory::Disulfide`. The bond category is set during
`CovalentTopology::Resolve` based on SG-SG covalent distance and
participating elements. Already in place; this slice consumes it.

The disulfide pair list passed to GenerateLeapScript uses 1-based PRMTOP
residue indices (LEaP convention). The mapping converts extractor
0-based residue indices into PRMTOP 1-based, accounting for caps that
shift the count.

### Determinism

Same Protein + same policy + same flat_table_path = same generated PDB,
same generated LEaP script, same PRMTOP (modulo tleap nondeterminism in
internal float arithmetic, which is reproducible on the same machine with
the same AmberTools build). The byte-equality test in step 4 is
machine-local; it confirms our generators are deterministic, not that
tleap is.

### Test taxonomy

```text
Pure-logic tests (no tleap, no AmberTools):
    AmberChargeResolver suite
    AmberLeapInput suite
    AmberPreparedChargeSource hard-precondition suite

Integration tests (require AMBERHOME / tleap):
    AmberPreparedChargeSource roundtrip suite
    Tagged AMBER_REQUIRES_TLEAP and SKIPped if RuntimeEnvironment::TleapPath()
    is empty.

Existing acceptance tests:
    ChargeFF14SB suite, ApbsFF14SB suite, PrmtopChargeTest, OrcaRunTest,
    FleetLoaderTest.ChargesReturned. All must remain green.

Blessed NPY parity:
    tests/golden/blessed/fleet/* unchanged for proteins on the existing
    flat-table or upstream-PRMTOP path.
```

## Open Decisions

### O1. Default `AmberPreparationPolicy` for `--pdb` / `--protonated-pdb`

```text
Conservative: FailOnUnsupportedTerminalVariants
    Behaves identically to current loud-fail. Capping is opt-in per run.

Permissive: UseCappedFragmentsForUnsupportedTerminalVariants
    Calibration runs that hit unsupported terminals "just work" with a
    documented capping decision in BuildContext.
```

Recommend conservative default. Permissive becomes the curated-corpus
default once the audit (see step 6 acceptance) confirms the cap count
and per-protein impact.

### O2. `--orca` strict-PRMTOP requirement — LOCKED 2026-04-29

`--orca` and `--mutant` are AMBER prmtop paths by design. The PRMTOP is
authoritative; nothing else is acceptable on those load paths.

Decision:

```text
BuildFromOrca rejects (returns hard error) when files.prmtop_path is
empty or the file does not exist. There is no fall-through to the
ff14SB flat table on the ORCA paths.
```

Mechanically the check already exists at OrcaRunLoader.cpp:304-314 and
must remain a hard fail when the resolver is wired in step 2. The
resolver's branch 1 (upstream PRMTOP) is a no-op for ORCA loads in the
sense that the loader has already guaranteed the precondition by the
time the resolver is called.

### O3. Cap geometry — NOT A QUESTION 2026-04-29

ff14SB is a fixed-charge force field. tleap reads atom names from the
input PDB, looks up CHARGE and RADII values from the residue library
template, and writes those values to the prmtop verbatim. Atom positions
do not enter the charge or radius computation. mbondi2 radii are also
library lookups, not position-dependent.

Therefore cap atom positions can be any physically reasonable values
that pass tleap's input validation. Ideal AMBER geometry from standard
bond lengths and angles is sufficient permanently; there is no future
"refinement" slice for this. The DFTs and any other geometry-sensitive
downstream step are not affected because cap atoms are dropped from the
result before charges leave the AmberPreparedChargeSource boundary.

### O4. NME as default C-terminal cap — LOCKED 2026-04-29

Decision: `UseCappedFragmentsForUnsupportedTerminalVariants` uses ACE
for N-terminal caps and NME for C-terminal caps.

Rationale:

```text
NME (-C(=O)-N(H)-CH3) mimics "a next ALA-like residue", so the original
residue keeps its INTERNAL backbone chemistry as if it were mid-chain.
This is the conservative choice: the cap exists only to obtain
INTERNAL charges from the stock library; we do not want the cap to
change the chemistry of the original residue.

NHE (-C(=O)-N(H)(H)) is a primary amide. It naturally occurs in some
peptides (amidated termini, e.g., neuropeptides), but for proteins
where the C-terminus is a real chain end with an unsupported variant,
using NHE would mis-represent the chemistry. Not used in this slice.
```

If a future input is a genuinely amidated peptide, that is a different
case and gets its own enum value on `AmberPreparationPolicy`. Out of
scope here.

## Dependencies and Required Environment

```text
AmberTools (tleap) installed and AMBERHOME exported, OR
RuntimeEnvironment::TleapPath() returns a valid path.
    On batcave: /home/jessica/micromamba/envs/mm/bin/tleap.
    Tests gated by AMBER_REQUIRES_TLEAP skip when absent.

Existing dependencies retained:
    cif++ (PDB parsing in PdbFileReader)
    libgromacs (TPR path, unchanged)
    HighFive, Eigen, etc. (unchanged)
```

## Rollback Strategy

Each step is its own commit. If step N fails NPY parity unexpectedly:

```text
git revert HEAD  (or git revert <N> if N is older)
```

Steps 1–3 cannot fail NPY parity by construction (no behavior change).
If they do, that indicates a bug in the refactor itself; revert and fix.

Steps 4–6 cannot affect NPY parity for proteins on the existing path
under FailOnUnsupportedTerminalVariants default. If they do, revert.

If step 6 capping produces incorrect charges for any test case: revert
step 6, leave 1–5 in. The flat-table loud-fail behavior remains.

## Post-Slice Cleanup Queue (NEXT, not "later forever")

The following cleanup is queued to run **immediately after** this AMBER
slice closes (step 6 acceptance gate green). It is named here so it does
not slip into "deferred forever."

### Sequencing (locked 2026-04-29)

```text
PHASE 1 — Build the COMPLETE chemistry+geometry substrate FIRST.
  N1.A through N1.G land in full, before any existing calculator is
  touched.

  "Complete" is the load-bearing word: every calculable typed field
  gemmi/CCD/RDKit/our-parsing can produce is on the topology, the
  projection surface is written, the MD-condition fields are typed,
  the pH/ionization layer is exhaustive. We squeeze out all the
  information this layer will ever hold before moving on.

  New calculators are allowed in this phase ONLY when they encode
  the new typology with geometry from day one (Phase 1.D's ring
  geometry, peptide ω, disulfide χss, rotamer classification).
  Those are one-offs that consume the new substrate natively, not
  migrations of existing calculators. Existing calculators are
  untouched in Phase 1.

  Why: a calculator migrated against an incomplete substrate will
  need to be migrated again when the substrate completes. Doing
  it twice doubles the review surface and doubles the chance of
  drift. The discipline is: substrate first, complete; calculators
  later, once.

  If we discover during Phase 2 that the substrate is missing a
  field a calculator needs, we extend Phase 1 — we never migrate
  a calculator against an incomplete substrate.

PHASE 2 — Calculator-by-calculator migration.
  ONLY after Phase 1 is closed. Each existing calculator migrates
  ONCE, to the final typed surface that already exists. Per-
  calculator drift report. Per-calculator review (drift OK,
  surprise not). Calculators that do not need richer information
  stay on the legacy surface indefinitely; no forced migration.

OpenBabel exit — CONSEQUENCE of Phase 2 progress.
  When the last existing calculator that depends on OpenBabel-derived
  data has migrated (or has been confirmed to remain on the legacy
  surface but with its OpenBabel-touching paths replaced by typed
  reads on the now-complete substrate), OpenBabel falls out as unused
  and the build dependency is removed. There is no "swap OpenBabel"
  phase; OpenBabel is retired by being unread.

PHASE 0 (sequencing-wise next, before Phase 1) — AMBER trajectory
  loader. Seals the "we are AMBER now" stance across every input
  path.

  GROMACS doctrine (locked 2026-04-30): all GROMACS interaction in
  this project is via the linked C++ API (libgromacs / xdrfile /
  read_tpx_state / etc.). Subprocess calls to `gmx dump`,
  `gmx trjconv`, `gmx editconf`, etc. are NOT acceptable — even
  when they would be "easier." The single existing subprocess
  GROMACS path (`GmxTprChargeSource` in src/ChargeSource.cpp,
  shells out to `gmx dump`) is legacy/quarantined and is being
  superseded by this slice. PHASE 0 must use C++ APIs for every
  GROMACS read; output writing (PDB/mmCIF frame extraction) goes
  via gemmi C++, not via shelled-out `gmx trjconv`.

  Why locked: the only legitimate "external call" reason is a
  physics need that the C++ API cannot satisfy (none currently
  identified). "Easier to wire" is not a physics reason and not
  acceptable.
  After PHASE 0 lands:
      --pdb                 → AMBER ff14SB ✓
      --protonated-pdb      → AMBER ff14SB ✓
      --orca / --mutant     → AMBER prmtop ✓
      --trajectory          → AMBER prmtop + AMBER frame stream  ← new
  All five input methods are AMBER-native. The CHARMM/GROMACS-from-
  TPR loader becomes quarantined legacy.

  Why before Phase 1: every Phase 1 substrate field can then have
  test coverage that spans all five input methods, exercising the
  same substrate built from each of them. This is the "we are AMBER
  now and the AMBER bit is wrapped" state.

  Test-fixture swap: the existing CHARMM trajectory test fixture
  (currently driving FleetLoaderTest cases) is replaced with an
  AMBER trajectory fixture as part of this phase. The FleetLoaderTest
  cases are reworked against the new fixture. The legacy CHARMM-from-
  TPR loader path stays in the codebase as quarantined fallback;
  it does not gate any new test.
```

### Crystal projection rule (load-bearing)

This rule is what prevents the typed-topology cleanup from collapsing
back into a pile of strings. **Read it as part of every N1 commit
review.**

```text
A projection in this design is a PURE FUNCTION on the typed substrate
fields. It is never a stored string. Every call recomputes the
projected name from the typed fields:

    std::string LegacyAmberTopology::IupacAtomName(size_t atom_index)
    const {
        const auto& f = atom_typed_fields_[atom_index];
        return IupacFormatter::Format(
            f.locant, f.branch_index, f.diastereotopic_index,
            f.methyl_methylene_class, f.polar_h_kind);
    }

NOT acceptable:

    std::string LegacyAmberTopology::IupacAtomName(size_t atom_index)
    const {
        return cached_iupac_names_[atom_index];   // antipattern
    }

The discipline:
    - Typed fields are the single source of truth.
    - Names are computed from them at call time.
    - A cached name is a derived value that can drift from the
      substrate; if the substrate is the truth, drift must be
      forbidden by construction.

If at any future review someone proposes "let's cache the IUPAC names
for performance," the answer is: profile first; if name generation is
hot, optimise the projection function (table lookup, etc.) but never
store the result. Drift from substrate is methodologically intolerable.

Output writers (H5 emitters, methods text generators, BMRB joiners)
call the projection at write time and consume the returned string
locally. Nothing stores it on the topology side.
```

### N1. Build the complete typed picture (LegacyAmberTopology + Conformation calculators)

**Goal.** Pre-compute and store, before any NMR calculator runs, every
calculable fact about the protein's chemistry and geometry that
calculators or downstream statistics could plausibly want. Calculators
read from this rich substrate; they do not re-derive any of it.

**Boundary discipline (load-bearing).** Where each fact lives is
governed by what it depends on:

```text
INVARIANT under our fixed-protonation MD
    → lives on LegacyAmberTopology (per-protein)
    → computed before FinalizeConstruction
    → never depends on atom positions

CONFORMATION-DEPENDENT (geometry; "stuff moves")
    → lives on ProteinConformation
    → computed by calculators (GeometryResult-style)
    → may depend on atom positions

A field on LegacyAmberTopology that depends on positions is a
category error. A field on ProteinConformation that depends only on
chemistry is also a category error. The split is enforced by where
the field lives.
```

**Layered phases:**

#### N1.A — Substrate parser pass (invariant; no library; LegacyAmberTopology)

Parse AMBER atom names into typed enums. This is the "names go in at
the boundary; typed fields come out" step that lets every other layer
read structured facts instead of strings.

```text
Per atom (typed enum or integer):
    locant                  (α, β, γ, δ, ε, ζ, η, ...)
    branch_index            (1, 2 for HD11/HD12 etc.)
    diastereotopic_index    (1, 2, 3 for methyl Hs / methylene Hs)
    methyl_methylene_class  (Methyl / Methylene / Methine / NotApplicable)
    polar_h_kind            (Amide / Hydroxyl / Sulfhydryl /
                             Charged / NotPolar)
```

No external library; pure C++ on the AMBER atom name string at the
construction boundary, then enum values forever after.

#### N1.B — CCD ingest (invariant; gemmi; LegacyAmberTopology)

```text
Per atom:
    stereo_config           (R / S / pro-R / pro-S / NotChiral)
    formal_charge           (chemistry, signed int)
    h_bond_donor / acceptor (chemistry-defined bool)
    electronegativity, atomic_polarizability, vdw_radius
                            (literature constants per element)

Per bond:
    bond_order              (Single / Double / Aromatic / Peptide)
    is_aromatic             (cross-validated with our ring set)
    is_in_ring              (CCD bond ring-membership flag)

Per residue:
    pKa per ionizable group (CCD or literature pKa table)
    functional_groups       (typed enum list with per-group atom set)
    ring_class verification (cross-checks our RingTypeIndex)
```

**Sourced from**: gemmi reading CCD (`chem_comp_atom`, `chem_comp_bond`,
`pdbx_chem_comp_descriptor`).

#### N1.C — pH / ionization layer (invariant for our MD; LegacyAmberTopology)

For our MD setup protonation does not change during a trajectory, so
pH-determined state is invariant once pH is supplied at construction.

```text
Per ionizable residue:
    protonation_state_at_pH       (already typed via
                                   protonation_variant_index)
    fraction_protonated_at_pH     (Henderson-Hasselbalch on residue pKa)
    charge_at_pH                  (formal charge under dominant state)

Per protein:
    net_charge_at_pH              (sum across residues; near-integer)
    isoelectric_point             (computed from pKa list)
```

`ProteinBuildContext::protonation_pH` already carries the pH input.
This layer reads it and populates the rest deterministically.

#### N1.D — Per-conformation geometric topology (calculator output;
                                                  ProteinConformation)

**These are calculators, not topology** — they depend on atom
positions, which move in MD. They follow the existing GeometryResult
pattern: read positions, compute typed values, store on
ProteinConformation alongside the existing per-frame fields
(`ring_geometries`, `bond_lengths`, `bond_directions`, `bond_midpoints`).

```text
Per ring (per frame):
    centroid (Vec3)
    normal vector (Vec3)
    radius (mean atomic distance from centroid)
    planarity_rms_deviation

Per ring system (fused; per frame):
    aggregate centroid, normal, atom-set membership

Per peptide bond (per frame):
    omega dihedral (planarity check)
    cis_trans_flag (typed enum)

Per disulfide (per frame):
    sg_sg_distance
    chi_ss dihedral

Per residue chi (per frame; partly already exists):
    chi values
    rotamer_classification (Gauche+/Gauche-/Trans, typed enum)
        — derived from chi value via a typed classifier; the Dunbrack
          rotamer LIBRARY data lives on LegacyAmberTopology as
          parameters
```

These land as new calculators (or extensions of `GeometryResult`).
The Dunbrack rotamer library data itself is invariant parameter data
and lives on `LegacyAmberTopology`; the per-frame classification is
the calculator output that uses the library + the chi value.

#### N1.E — RDKit perception (invariant chemistry; LegacyAmberTopology;
                                optional / heavier)

```text
Per atom:
    hybridisation                 (cross-validate CCD; sp/sp2/sp3 enum)
    is_aromatic                   (cross-validate CCD)
    gasteiger_partial_charge      (theoretical, optional)
    mmff94_partial_charge         (theoretical, optional)

Per residue / chain:
    canonical_smiles              (output projection helper, optional)
    functional_groups via SMARTS  (cross-validate N1.B)

Per bond:
    bond_dipole_moment, bond_polarizability
                                  (literature lookup keyed on bond_order
                                   + endpoint elements)
```

RDKit perception runs at construction; the resulting typed values are
stored on LegacyAmberTopology. Theoretical partial charges are NEVER
substituted for `ForceFieldChargeTable::PartialChargeAt` — they're
parallel data for cross-validation, ML features, etc.

This phase is **optional / heavier** because RDKit is a substantial
build dependency. We add it only if downstream value justifies the
weight.

#### N1.F — Projection surface (functions on LegacyAmberTopology)

Pure functions, never cached, computed every call from typed fields:

```text
amber_atom_name(atom_index)         : string
iupac_atom_name(atom_index)         : string
bmrb_atom_name(atom_index)          : string
amber_residue_name(residue_index)   : string
iupac_residue_name(residue_index)   : string
canonical_smiles(chain_index)       : string  (RDKit, when enabled)
```

Each is a thin formatter on the typed substrate. See "Crystal
projection rule" above. Output writers call these at write time;
nothing on the topology side stores the result.

#### N1.G — The substrate as MD-condition specification

The typed substrate that N1.A through N1.F builds out is **exactly
the protein-condition specification an MD setup consumes**.

```text
Protein-condition fields covered by the substrate:
    pH                                    (already on BuildContext)
    per-residue protonation_state_at_pH   (N1.C)
    net_charge_at_pH                      (N1.C)
    isoelectric_point                     (N1.C)
    force_field_name                      (already on BuildContext)
    temperature                           (TO ADD as a typed field on
                                           ProteinBuildContext alongside
                                           protonation_pH)
    pKa_source                            (CCD or curated table; cite)
    counter_ion_neutralization_rule       (derivable from net charge
                                           when a rule is documented)

Run-side fields that DO NOT belong on the substrate (trajectory
metadata, not protein):
    water model (TIP3P / TIP4P-Ew / OPC / SPC/E)
    box dimensions, PBC
    equilibration protocol (stages, durations)
    integrator and timestep
    thermostat type and coupling time
    pressure-coupling type
    frame interval
    random seed
```

**Two methodological consequences:**

```text
1. Round-trip integrity.
   Loading a trajectory back through the cleanup pipeline regenerates
   the same substrate identically. No pH guessing, no protonation
   re-derivation, no force-field-choice ambiguity. The trajectory's
   protein side is fully reconstructible from typed state — the
   round-trip discipline thesis reproducibility wants.

2. Cross-trajectory comparability becomes typed.
   When two trajectories ran at different conditions, the substrate of
   each carries the difference as typed values. Comparisons and plots
   read the substrate directly; nobody parses run notes to know which
   is which.
```

**Implementation cost.** Small. Add one field
(`ProteinBuildContext::temperature_kelvin`) and a citation slot for the
pKa source. The conceptual recognition is the load-bearing part — the
substrate IS the MD-condition spec, not coincidentally but by design.

### N2. Calculator-by-calculator migration

After N1 lands, calculators that benefit from the rich typed fields can
migrate at their own pace. Each calculator review either:

```text
- keeps reading from the legacy surface (Bond::category,
  Atom::bond_indices, Atom::pdb_atom_name) — fine, no change required
- migrates a specific read to the typed accessor (e.g. methyl-H
  averaging via is_methyl_H typed flag instead of pdb_atom_name
  string match) — produces a per-calculator drift report
```

Per-calculator drift reports follow the rule from earlier in this
document: **drift is OK; surprise is not**. Methods text cites the
migration date and the documented cause of any numerical change.

Calculators that do NOT need richer information stay on the legacy
surface indefinitely. There is no forced migration.

### N3. OpenBabel exit (consequence)

OpenBabel is retired when no remaining calculator (and no remaining
construction-time path) reads its outputs. This is a build-dependency
removal, not a feature change. By the time we get here, the typed
topology already provides everything calculators need; OpenBabel just
stops being linked.

### N4. Trajectory loading consolidation

**Status.** Queued for immediately after step 6 acceptance gate. Needs
**design discussion before implementation** — see "Completeness
constraint" below. This is not a one-session swap.

**Context.** `CovalentTopology::Resolve` currently uses OpenBabel
(`HAS_OPENBABEL=1` is always defined) to detect which atom pairs are
bonded; classification (Disulfide / PeptideCO / PeptideCN / Aromatic /
BackboneOther / SidechainOther) is then typed. For curated
AMBER-typed proteins — which is the project's input domain — every
intra-residue bond is determined by `AminoAcidType` templates, every
peptide bond is determined by sequential residue order within a chain,
and every disulfide is determined by inter-CYS SG-SG distance.

**Methodological note (defensible-by-citation).** AmberTools tutorials
prescribe disulfide detection via inter-CYS SG-SG distance. tleap reads
intra-residue bonds from its residue libraries (`amino12.lib`,
`aminont12.lib`, `aminoct12.lib`) — it does not perform geometric bond
perception. A typed bond resolver mirrors what AMBER itself does.

### Completeness constraint (the hard part)

Bond information is **not just an AMBER input** — it is **calculator
and statistics input**. `BondCategory`, `Atom::bond_indices`,
`Atom::parent_atom_index`, and the typed bond list feed:

```text
- BiotSavart and Haigh-Mallion bonded-weight terms.
- McConnell ring-current calculations (via aromatic bond perception).
- Bond anisotropy contributions.
- Peptide-bond-specific terms (PeptideCO, PeptideCN distinct).
- Mutation delta (matches across WT/mutant via bond context).
- Calibration ridge regression (input features depend on bond
  classification).
```

### Drift policy (corrected 2026-04-29)

**NPY drift is acceptable per phase, provided we characterise and
defend every difference.** Strict bit-identity is engineering hygiene;
*understanding* the drift is the science. The two correspond when
nothing changes, but the load-bearing principle is the second one.

```text
Drift is OK.
Surprise is not.
Every numerical change must have a named, defensible reason that
survives review.
```

The completeness constraint survives but is reframed: we do not want to
*silently lose* a bond OpenBabel was finding that calculators relied
on. But "completeness" no longer requires "identical numbers." It
requires:

```text
1. Every bond list difference is enumerated.
2. Every calculator output difference is traced to a bond list
   difference (or another typed-source difference).
3. Every difference has a documented reason — typed methodology
   citation, OpenBabel quirk, geometric edge case — that the methods
   text can defend.
```

### Per-phase drift reporting

Every cleanup phase produces a small drift report on the curated
corpus:

```text
Phase X drift report (curated corpus, N proteins):
  bond list differences:
    K bonds added vs OpenBabel — listed with per-protein incidence
    M bonds removed vs OpenBabel — listed with per-protein incidence
    each with cited reason ("typed sequential C-N peptide bond confirmed
      by distance check; OpenBabel had failed to detect at chain-break
      gap" / "OpenBabel quirk: spurious vdW-near pair classified as
      Single bond; typed resolver omits correctly" / etc.)
  calculator output drift:
    max |Δ| per kernel
    root-causing bond list change for each delta exceeding {threshold}
  decision: bless / investigate further / revert this phase
```

This drift report becomes part of the thesis-methodology trail.
Methods text cites it: "Bond classification was migrated from
chemistry-perception (OpenBabel) to typed templates (AmberTools-aligned
methodology). The migration produced K differences across the
curated corpus, all traceable to {documented cause classes}."

### Audit deliverables (concrete, pre-implementation)

Three deliverables before any cleanup-phase code lands:

```text
1. Field-consumer audit:
   Grep every read of Bond::is_rotatable, Bond::order,
   Bond::is_aromatic, Atom::bond_indices, Atom::parent_atom_index
   across src/. List which calculators consume which fields. This is
   the surface that must continue to populate.

2. Typed-source mapping:
   For every consumed field, name its typed source post-cleanup
   (AminoAcidType template, sequential C-N walk, SG-SG distance,
   ring lookup, chi-angle definition, etc.). If no typed source
   exists today, decide: (a) extend AminoAcidType / sibling table
   to cover it, (b) add a named geometric check, or (c) retire
   the field with a documented reason.

3. Surface-preservation classification:
   Per field, classify what's expected:
       SAME — typed source produces identical values for the
              curated corpus.
       DRIFT-DOCUMENTED — typed source produces different values
                         for documented reasons; methods text cites.
       NEW — typed source provides new information not previously
             available (richer typed fields per N2).
```

These three deliverables together let us answer "is the calculators-
unchanged promise realistic for this field?" honestly per field.

### Doctrine

```text
Drift with understanding is the bar.
Methods text cites the drift report.
A passing-but-undefended diff is a fail.
```

**Why this is a discussion before code, not "rip and replace."** The
audit must be done first because the calculator-by-calculator
migration ordering depends on which fields drift and how. We sequence
the easy-to-defend phases first; the harder ones get more design time.

### Candidate replacements (for the design discussion)

User-surfaced 2026-04-29: the candidates worth comparing against the
typed-only resolver are:

```text
RDKit         — open source, very mature bond perception, used widely in
                pharma / drug-discovery codebases. C++ API. Free. Heavy
                dependency. Best general perception quality outside
                commercial tools.

OpenFF        — Open Force Field Initiative; modern force-field-aware
                tooling with SMIRKS-based perception. Python-first, less
                mature C++ story. Free. Specifically designed for
                force-field workflows; aligns well with our pipeline
                framing but adds a new toolchain.

OpenEye       — commercial license, considered gold standard for bond
                / chemistry perception in pharma. Methods-defensible by
                citation. Cost barrier and licensing complications.

GAFF / antechamber  — already shipped in AmberTools (which we depend on
                via tleap). The AMBER-blessed route for chemistry outside
                stock residue templates. Same methodology citation as the
                rest of our AMBER pipeline. No new external dependency.
                Specifically designed for force-field-grade bond /
                charge perception of organic molecules.
```

Read-through: for the project's input domain (curated standard proteins
with stock residues), all four give the same answer on the protein
core. They diverge on edge cases (non-standard residues, ligands,
unusual crosslinks) — and our resolver verdict already loud-fails on
those today, so divergence on the edge cases doesn't load-bear here.

**Strong default candidate**: **GAFF / antechamber via AmberTools**.
Rationale:
```text
- Already a dependency (we use tleap).
- Same toolchain, same citation, same methodology shape as the rest of
  the AMBER pipeline.
- "Bonds detected through AmberTools antechamber + standard residue
  templates" is one defensible methods sentence end-to-end.
- No new build dependency.
```

**Alternate candidate**: **RDKit**, if the project wants a chemistry
perception layer independent of AmberTools. Rationale:
```text
- Most general-purpose; highest code quality outside commercial tools.
- Useful beyond bond detection (substructure matching, conformer
  generation, etc.) if other slices want it.
- New build dependency (~hundreds of MB).
```

**Out of scope for this discussion**: OpenEye (license/cost
incompatibility with the project's open-toolchain stance unless the user
specifically opens that door); OpenFF (immature C++ story; would force a
broader Python/C++ split for this layer alone).

The audit (sub-questions 1-3) plus the candidate comparison is what
the post-slice design session is for. No implementation choice yet.

---

## Known follow-up code items (review-flagged, deferred)

These were surfaced during the 2026-04-29 evening doc review. They are
real but non-urgent; tracked here so they're not lost.

```text
F1. Shell-quoting in tleap subprocess invocation.
    Location: src/AmberPreparedChargeSource.cpp around line 336.
    Today the command is built via string concatenation:
        cmd = tleap_bin + " -f " + script_path
              + " > " + log_path + " 2>&1"
    and passed to std::system. In current use the path components
    come from RuntimeEnvironment::TempFilePath
    (/tmp/nmr_shielding/<guid>_..., no spaces / no shell metachars),
    so this is robust by inputs. Defensive cleanup: route through
    posix_spawn or a small subprocess helper that takes argv and
    log redirection explicitly. Low priority; not a security issue
    in this codebase's threat model.

F2. Capping match ignores residue insertion code.
    Location: src/AmberLeapInput.cpp around line 178-182.
    BuildCappingPlan locates the failure-residue by chain_id +
    sequence_number; it does not check insertion_code. PDB inputs
    with ICODE'd residues sharing a sequence_number (100 / 100A /
    100B) are technically ambiguous — the first match wins.
    Impact: bounded — (unsupported terminal variant) AND
    (ICODE'd residue) is a rare combination. Real bug regardless;
    extend the match predicate to include insertion_code when this
    is touched next.

F3. PRMTOP parser duplication.
    AmberPreparedChargeSource.cpp duplicates the prmtop section
    readers (ReadPrmtopFlagLines / ReadPrmtopDoubles /
    ReadPrmtopInts / ReadPrmtopStrings) that also live in
    PrmtopChargeSource.cpp and OrcaRunLoader.cpp. Three copies of
    near-identical code. Unify into a shared src/PrmtopParser.{h,cpp}
    when convenient; not a correctness issue.
```

---

## Doctrine

```text
Source path lives in typed fields, not loader switches.
Resolver is the only construction point for AMBER ChargeSources.
tleap runs only when no upstream artifact exists AND the flat table
    cannot represent the protein.
Re-tleap is impossible by construction; an assert protects the rule.
Atom mapping by typed (residue_index, atom_name); cap atoms dropped by
    residue name, never by guess.
Generated PDB and LEaP script are deterministic functions of typed
    Protein state.
NPY parity is the regression gate at every step.
Provenance is typed and lives in BuildContext, available to methods text.
```
