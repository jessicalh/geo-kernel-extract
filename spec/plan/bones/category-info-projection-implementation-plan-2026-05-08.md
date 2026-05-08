# CategoryInfoProjection — implementation plan (2026-05-08)

**Status:** plan locked, zero implementation. Recovered and reframed
across the 2026-05-08 session after a GNOME UI crash. This is the
single durable record. If another crash hits before this lands, this
is the doc to read first.

This slice produces **the comprehensive per-atom categorical record**
— every invariant fact we know about each atom, in one structured
NPY, in one place. It exists to unblock codex's Stage 1 stats redo
and to replace the geometry-based atom matchup in
`MutationDeltaResult` with typed-identity matching. Downstream
calibration analysis pivots every stratification on these fields.

## 1. What landed pre-crash

### Two commits on 2026-05-08

| SHA | Time | Title | Notes |
|---|---|---|---|
| `6c4d13e` | 10:24 BST | BMRB nomenclature: add atom_nom.tbl reference + ignore molprobity_runs/ | Persists `references/bmrb_data/atom_nom.tbl` (BMRB hydrogen atom-naming cross-walk, IUPAC/IUB 1970 + Markley 1998). Bibliography section explains scope. `.gitignore` adds `molprobity_runs/`. |
| `85c1637` | 10:26 BST | FramePdbEmitter production validation: HDF5 + CRYST1 fixes + docs + MolProbity wrapper | HDF5 1.14↔1.10 mismatch fix, one-char CRYST1 transpose fix, USE_CASES doc, `tools/molprobity_validate.py`. All production-validated against the 1P9J 751-frame extract. |

### Two memory entries saved (in `-shared-2026Thesis/memory/`)

- `feedback_naming_input_output_asymmetry.md` — load-bearing
  architectural rule: input naming feeds physics (no tolerance,
  fail-loud); output projection feeds ML matching (logged fallbacks
  acceptable). NEVER glue them.
- `feedback_bmrb_chemistry_not_raw.md` — don't diff `chemistry.json`
  against naive BMRB `_Sample_condition_variable` reads.

### Side-products on disk (gitignored, regenerable)

- `molprobity_runs/1P9J_5801/` — 1P9J 751-frame PDB tree + MolProbity validation.
- `/shared/2026Thesis/paper_archive/1P9J_5801_15ns_optB_20260501/` — paper-archival 1P9J extract.

### Configuration in place

- `~/.nmr_tools.toml` has `bmrb_atom_nom = "/shared/2026Thesis/nmr-shielding/references/bmrb_data/atom_nom.tbl"`.
- `references/bmrb_data/atom_nom.tbl` committed, ~12 KB, ~400 entries.

## 2. Architectural rule

From `feedback_naming_input_output_asymmetry`:

> Input-side and output-side naming systems are NEVER glued together.

- **Input side** (loader → `Atom.pdb_atom_name` → `NamingApplicator`
  → typed substrate composition): wrong name corrupts physics
  silently. Fail-loud on unknown.
- **Output side** (typed substrate → IUPAC / BMRB / categorical
  record → NPY): errors cost data points, not physics.
  Logged-fallback acceptable.

`CategoryInfoProjection` lives in its own TU. Not a method on
`LegacyAmberTopology` (`OBJECT_MODEL.md:420` is explicit about this:
"It does NOT live on `Protein`, on `LegacyAmberTopology`, or on
`Atom`"). Output-only invocation, NPY-emission boundary, gated on
`HasAtomSemantic()`.

## 3. Locked decisions

- **Name: `CategoryInfoProjection`.** "Naming projection" undersold
  it — this is the comprehensive categorical identity record.
- **Single structured NPY per protein.** ~28 fields, mixed `S8` /
  `S4` / `S1` / `i1` / `i4`. ~50 bytes per atom × ~4000 atoms = ~200
  KB per protein. **One** `ArraySpec` entry, **one** SDK wrapper class.
- **One-shot emission per protein.** Not per-conformation; not
  per-frame. `WriteAllFeatures` is unchanged. Two new call sites in
  the entry points (non-trajectory entry + `Trajectory::Run` before
  the per-frame loop).
- **FramePdbEmitter shape.** Static methods, file-local state in
  `CategoryInfoProjection.cpp` anonymous namespace, deleted ctor + copy
  + assign. No instance, no Session field, no global free function.
  Configured once at startup from `Session::LoadFromToml`.
- **`bmrb_atom_name` and `iupac_atom_name` as separate columns**
  even though they're identical for the standard 20. Schema
  future-proofing: divergence (non-standard residues, pseudoatom
  Q-labels) becomes a row-level change, not a schema change.
- **Provenance enum: `Match / MissLogged` only.** No `DriftLogged`.
  Substrate fields are invariants — substrate-vs-`atom_nom.tbl`
  disagreement is a chemistry error, not a routine fallback. If it
  surfaces in real data, treat as a substrate bug.
- **`S8` strict-length fail-loud guard at emission.** Any
  `atom_name` or `residue_3letter` exceeding 7 chars aborts with a
  diagnostic naming the offending row. Per-PATTERNS §9: `fprintf +
  std::abort`, no exceptions.
- **MutationDeltaResult matchup rewrite is in scope.** Replace
  greedy-spatial+element matching at `MutationDeltaResult.cpp:278-330`
  with typed-identity matching using
  `(residue_index, AtomMechanicalIdentity)`. Both pieces (the NPY +
  the rewrite) gate codex's Stage 1 stats redo on 723 proteins.
- **Per-residue fields denormalized to per-atom.** Residue 3-letters,
  variant index, terminal state copied to every atom in the
  structured record for join-free filtering.
- **Main-thread implementation, no in-codebase agents.** Per the
  FramePdbEmitter precedent.

## 4. File plan

### Add

| Path | Purpose |
|---|---|
| `src/CategoryInfoProjection.{h,cpp}` | FramePdbEmitter-shaped projection class. Parses `atom_nom.tbl` once at `Configure`, builds the per-atom record on `WriteFeatures`, emits one structured NPY. Holds static maps for canonical 3-letter / 1-letter projections. Owns the miss-log diagnostic. |
| `tests/test_category_info_projection.cpp` | Construction (parses atom_nom.tbl), per-atom record fields against 1ubq / 1P9J / Trp-containing fixture, variant 3-letter mapping, prochiral stereo cases including Gly inversion, fallback path, S8 length-guard fail-loud. |
| `tests/test_mutation_delta_typed_match.cpp` (or extend existing) | The typed-identity match correctness test: a known WT/mutant pair where greedy-spatial misbinds, verify typed match binds correctly. |

### Modify

| Path | Change |
|---|---|
| `src/RuntimeEnvironment.{h,cpp}` | Add `BmrbAtomNom()` static accessor + `bmrb_atom_nom_` field + TOML load. Same pattern as `Mopac()`, `Tleap()`, `Ff14sbParams()`. |
| `src/Session.cpp` | After `RuntimeEnvironment::Load()`: `CategoryInfoProjection::Configure({.atom_nom_tbl = RuntimeEnvironment::BmrbAtomNom()})`. One line. No new Session field. |
| `src/nmr_extract.cpp` (or whichever entry point handles non-trajectory paths) | Call `CategoryInfoProjection::WriteFeatures(protein, output_dir)` once after Protein construction, before the per-conformation work. |
| `src/Trajectory.cpp` | Same call once before the per-frame loop in `Trajectory::Run`. |
| `src/MutationDeltaResult.cpp` | Replace lines 278-330 (greedy-spatial+element atom matching) with typed-identity matching via `LegacyAmber().IdentityAt(ai)`. Bijection automatic; mismatches surface fail-loud. |
| `CMakeLists.txt` | Add `src/CategoryInfoProjection.cpp` to `nmr_shielding` source list. |
| `python/nmr_extract/_catalog.py` | One new `ArraySpec` entry for the structured NPY. The wrapper SDK exposes typed Python enums for the int8 columns and `bytes`→`str` decode for the S8 columns. |
| `python/nmr_extract/Protein.py` (or new `_category_info.py`) | New SDK class wrapping the structured NPY. Named accessors with typed enum returns. |

### Order

1. `RuntimeEnvironment::BmrbAtomNom()` accessor.
2. `CategoryInfoProjection.{h,cpp}` — class + atom_nom.tbl parser + canonical 3-letter / 1-letter maps + record builder + structured-NPY writer + tests.
3. Two one-shot call sites in entry points.
4. `MutationDeltaResult` matchup rewrite + test.
5. SDK catalog + Python wrapper + Python tests.
6. `Session::LoadFromToml` calls `Configure`.
7. Build + ctest + commit (single bundled commit).

## 5. CategoryInfoProjection API

```cpp
namespace nmr {

class Protein;
enum class AminoAcid : uint8_t;

class CategoryInfoProjection {
public:
    struct Config {
        std::filesystem::path atom_nom_tbl;  // empty = inert (no atom_nom.tbl-driven names)
    };

    static void Configure(Config config);              // called once at startup
    static int  WriteFeatures(const Protein& protein,
                              const std::string& output_dir);
    static void Reset();                                // for tests
    static bool IsActive();

    // Per-atom queries — for tests and downstream callers (viewer / h5-reader).
    // All const; all read-only on the protein/topology.
    static std::string IupacAtomName(const Protein&, std::size_t atom_index);
    static std::string BmrbAtomName(const Protein&, std::size_t atom_index);
    static std::string AmberResidueThreeLetter(AminoAcid type);   // CYX / HID / etc.
    static std::string IupacResidueThreeLetter(AminoAcid type);   // CYS / HIS / etc.
    static char        ResidueOneLetter(AminoAcid type);

    // Diagnostic — total miss tally for this projection's lifetime.
    static const std::map<std::string, int>& MissLog();

    CategoryInfoProjection() = delete;
    CategoryInfoProjection(const CategoryInfoProjection&) = delete;
    CategoryInfoProjection& operator=(const CategoryInfoProjection&) = delete;
};

}  // namespace nmr
```

Internal state (parsed `atom_nom.tbl`, miss log, configured-flag) lives
in an anonymous namespace inside `CategoryInfoProjection.cpp`. Same
idiom as `FramePdbEmitter`.

## 6. The structured NPY schema

Filename: `atoms_category_info.npy` (or shorter — open).

NumPy structured dtype, shape `(N,)` where N = atom count:

```python
np.dtype([
    # ── Identity ──────────────────────────────────────────────────
    ('atom_index',                 'i4'),   # 0..N-1; row index, here for explicit join
    ('residue_index',              'i4'),   # 0..R-1
    ('element',                    'i1'),   # Element enum

    # ── Atom names across naming systems ──────────────────────────
    ('amber_atom_name',            'S8'),   # canonical AMBER ff14SB (= post-NamingApplicator pdb_atom_name)
    ('iupac_atom_name',            'S8'),   # IUPAC standard from atom_nom.tbl
    ('bmrb_atom_name',             'S8'),   # BMRB (== IUPAC for standard 20; placeholder for divergence)

    # ── Per-residue fields, denormalized to per-atom ──────────────
    ('amber_residue_3letter',      'S4'),   # CYX / HID / HIE / ASH / GLH / LYN / ARN / TYM / canonical
    ('iupac_residue_3letter',      'S4'),   # canonical (CYS / HIS / ASP / GLU / LYS / ARG / TYR / standard)
    ('bmrb_residue_3letter',       'S4'),   # canonical (== IUPAC)
    ('residue_1letter',            'S1'),
    ('residue_type',               'i1'),   # AminoAcid enum (the variant-specific one)
    ('residue_variant_index',      'i1'),   # AMBER protonation variant index
    ('terminal_state',             'i1'),   # Internal / NtermCharged / NtermNeutral / CtermDeprotonated / CtermProtonated

    # ── Mechanical identity (the AtomMechanicalIdentity tuple) ───
    ('locant',                     'i1'),   # Greek letter
    ('branch_outer',               'i1'),   # CG1 vs CG2; CD1 vs CD2; etc.
    ('branch_inner',               'i1'),   # second-level branching on H atoms
    ('di_index',                   'i1'),   # Position2 / Position3 on prochiral methylenes
    ('backbone_role',              'i1'),

    # ── Chemistry / topology classification ──────────────────────
    ('prochiral',                  'i1'),
    ('planar_group',               'i1'),
    ('planar_stereo',              'i1'),
    ('polar_h_kind',               'i1'),
    ('ring_position_primary',      'i1'),
    ('ring_position_secondary',    'i1'),   # second ring for Trp bridgeheads
    ('pseudoatom_kind',            'i1'),   # None / M / Q / R
    ('in_super_group',             'i1'),   # bool: also in a QG/QD/QH/QR aggregator
    ('aromatic',                   'i1'),
    ('formal_charge',              'i1'),
    ('is_exchangeable',            'i1'),

    # ── Provenance for any external table lookup ──────────────────
    ('iupac_naming_provenance',    'i1'),   # Match / MissLogged
    ('bmrb_naming_provenance',     'i1'),   # Match / MissLogged
])
```

### Where each field comes from

- `atom_index`, `residue_index`, `element` — direct from `Atom`.
- `amber_atom_name` — `Atom.pdb_atom_name` (already canonicalized by `NamingApplicator` to AMBER ff14SB).
- `iupac_atom_name`, `bmrb_atom_name`, `iupac_naming_provenance`, `bmrb_naming_provenance` — `atom_nom.tbl` lookup, keyed by `(residue_1letter, atom_name)`.
- `amber_residue_3letter` — small switch on `AminoAcid` enum; variants stay as CYX/HID/HIE/HIP/ASH/GLH/LYN/ARN/TYM.
- `iupac_residue_3letter`, `bmrb_residue_3letter` — same switch but variant residues collapse to canonical (HID/HIE/HIP → HIS, etc.).
- `residue_1letter` — small switch on `AminoAcid` enum.
- `residue_type` — `Residue.type` cast to int8.
- `residue_variant_index` — `Residue.protonation_variant_index`.
- `terminal_state` — `Residue.terminal_state`.
- `locant`, `branch_outer`, `branch_inner`, `di_index`, `backbone_role` — `protein.LegacyAmber().SemanticAt(ai).{...}` (the `AtomMechanicalIdentity` tuple).
- `prochiral`, `planar_group`, `planar_stereo`, `polar_h_kind`, `ring_position_primary`, `ring_position_secondary`, `aromatic`, `formal_charge`, `is_exchangeable` — direct copy from `AtomSemanticTable`.
- `pseudoatom_kind` — `SemanticAt(ai).pseudoatom.kind`.
- `in_super_group` — `SemanticAt(ai).pseudoatom.in_super_group` cast to int8.

## 7. Codex Stage 1 stratifications enabled

Every stratification becomes a direct numpy boolean mask, no joins:

```python
# per-atom-type
hb_atoms        = (atoms.element == ElementId.H) & (atoms.locant == Locant.Beta) & (atoms.di_index != DiIndex.None)
backbone_n      = atoms.backbone_role == BackboneRole.Nitrogen
sidechain_n     = (atoms.element == ElementId.N) & (atoms.backbone_role == BackboneRole.None)

# per-residue
phe_aromatic_c  = (atoms.iupac_residue_3letter == b'PHE ') & (atoms.aromatic == 1) & (atoms.element == ElementId.C)
his_imidazole_n = (atoms.residue_type == AminoAcid.HID) & (atoms.polar_h_kind == PolarHKind.ImidazoleNH)
ash_carboxyl_h  = (atoms.residue_type == AminoAcid.ASH) & (atoms.polar_h_kind == PolarHKind.CarboxylOH)

# pseudoatom-based aggregation
all_methyl_h    = atoms.pseudoatom_kind == PseudoatomKind.M
super_aggregate = atoms.in_super_group == 1   # Val QG, Leu QD, Arg QH, Phe/Tyr QR

# ring-position chemistry
trp_bridgeheads = atoms.ring_position_secondary != RingPos.None
phe_ipso        = (atoms.ring_position_primary == RingPos.Ipso) & (atoms.iupac_residue_3letter == b'PHE ')

# variant-specific
cyx_disulfide_partners = atoms.residue_type == AminoAcid.CYX
lyn_amine_h            = atoms.residue_type == AminoAcid.LYN
```

## 8. MutationDeltaResult matchup rewrite

### What today's matchup does (`MutationDeltaResult.cpp:278-330`)

Greedy-nearest-neighbor + element filter + bijection:

```cpp
struct Claim { size_t wt_index = SIZE_MAX; double distance = 1e30; };
// for each WT atom:
//   find mutant atoms within radius, same element only
//   greedy: claim closest, kick out previous claimer if newcomer is closer
```

No residue index. No atom name. No role. Position + element. A
rotated sidechain Hβ can quietly bind WT's Hα tensor against mut's Hβ
tensor; rings around the mutation site move and atoms within them
rebind incorrectly. The DeltaSummary that feeds calibration sees
these as legitimate per-element bins. This is the "cheesy" matchup.

### Replacement

Match by `(residue_index, AtomMechanicalIdentity)`:

```cpp
for (size_t wi = 0; wi < wt_count; ++wi) {
    if (wt_residue_is_mutation_site(wi)) continue;     // mutation site is intentionally unmatched
    const auto& wt_id = wt_protein.LegacyAmber().IdentityAt(wi);
    size_t wt_ri = wt_protein.AtomAt(wi).residue_index;
    size_t mi = mut_protein.FindAtomByIdentity(wt_ri, wt_id);
    if (mi == SIZE_MAX) {
        // Surface as named miss: residue/identity that exists in WT but not mutant.
        // Fail-loud: this should not happen for non-mutation residues.
    }
    // ... compute delta for matched pair
}
```

`Protein::FindAtomByIdentity(residue_index, AtomMechanicalIdentity)`
is the helper to add — walks the residue's atoms and returns the
unique atom whose substrate identity matches. Bijection automatic
because the identity tuple is unique per residue (verified by the
existing tests/topology tables work).

### What this enables

- Codex's Stage 1 stats redo gets matched-tensor pairs that are
  actually corresponding atoms.
- Geometric drift (rotamer flips between WT and mutant) doesn't
  contaminate per-element / per-atom-type bins.
- Mismatched cases surface as named errors instead of silent
  rebinds. Real chemistry mismatches (e.g. mutation alters bonding
  unexpectedly) become visible.

## 9. Test plan

`tests/test_category_info_projection.cpp`:

1. **Construction.** `Configure({.atom_nom_tbl = ".../atom_nom.tbl"})` parses without error.
2. **1ubq round-trip.** Standard 20 residues, basic backbone + sidechain naming.
3. **1P9J coverage.** Exercises HID/HIE/HIP variants, ASH/GLH if present, terminal caps, prochiral methylenes outside Gly.
4. **Trp-containing fixture.** Exercises `ring_position_secondary` (CD2/CE2 bridgeheads in both 5-ring and 6-ring).
5. **Variant 3-letter mapping.** HID/HIE/HIP → HIS, CYX/CYM → CYS, ASH → ASP, GLH → GLU, LYN → LYS, ARN → ARG, TYM → TYR.
6. **Gly Hα inversion.** HA2 = pro-R, HA3 = pro-S (the known inversion).
7. **Fallback / miss path.** A non-standard atom name returns empty string; `MissLog()` records the miss with appropriate provenance value.
8. **S8 length guard.** A synthetic 8+ char atom_name aborts at emission with a diagnostic naming the offending row.

`tests/test_mutation_delta_typed_match.cpp` (or extend
`tests/test_mutation_delta.cpp` if it exists):

1. **Bijection.** On a known WT/mutant pair, every non-mutation atom matches exactly one mutant atom by typed identity.
2. **Mutation site exclusion.** Atoms at the mutation site are intentionally unmatched.
3. **Robustness vs. rotamer drift.** Construct a fixture where rotamer flip would have caused the greedy matcher to misbind; verify typed match binds correctly.
4. **Fail-loud on missing identity.** Synthetic mutant missing an expected atom (e.g. CB removed) surfaces as a named error.

## 10. Out of scope (deferred)

- **Pseudoatom Q-label string projection.** BMRB's
  `pseudoatom_lib.tbl` is a separate cross-walk. The substrate's
  `pseudoatom_kind` + `in_super_group` boolean + the per-residue +
  locant context lets codex derive QG/QD/QH/QR letters in Python if
  needed; baking them as a separate column is its own slice.
- **AMBER-side projection cross-walk.** The AMBER name is already
  the substrate's authoritative form (`Atom.pdb_atom_name`). No
  separate projection module needed.
- **H5 emission.** `TrajectoryProtein::WriteH5` does not get the new
  fields in this slice. The structured NPY is the canonical store
  for now; an H5 parallel can be a follow-up if Stage 2 calibration
  analysis needs it.
- **Non-standard residues.** atom_nom.tbl covers the standard 20.
  Modified residues, ligands, cofactors → emit `MissLogged` with the
  raw `pdb_atom_name` as fallback in IUPAC/BMRB columns.

## 11. Constraints to respect

- **`§Adapter, wrapper, or bridge classes`** (PATTERNS) —
  `CategoryInfoProjection` is none of those. It's a named projection
  operation, sized appropriately. The per-atom queries are convenience
  for tests; the production hot path is `WriteFeatures` writing a
  structured NPY.
- **`§Naming boundary: cross once, then typed objects forever`**
  (PATTERNS) — the projection IS at the output boundary. Strings live
  inside `CategoryInfoProjection.cpp` and `atom_nom.tbl`; they leave
  only as `S8` bytes in the NPY. Calculators continue to read the
  typed substrate.
- **`§Return codes, not exceptions`** (PATTERNS §9) — S8 length-guard
  violation aborts with `fprintf(stderr, "FATAL: ...") + std::abort()`.
  No exception hierarchy, no `std::expected`.
- **`§The utility namespace`** (PATTERNS) — no `Util.h`,
  `Helpers.cpp`, or `util::` namespace. The projection is its own
  named concept with a class.
- **`§Singleton guarantee`** (PATTERNS §3) — does not apply here;
  the project's "singleton" concept is per-Result-per-Conformation.
  `CategoryInfoProjection` uses static-class file-local state, the
  `RuntimeEnvironment` / `FramePdbEmitter` pattern for process-wide
  utilities.
- **`feedback_no_attach_lifecycle_for_invariant_data`** — the
  projection holds invariant configuration (parsed atom_nom.tbl
  table) in static storage; no `Attach` / `Has` / `Optional`
  lifecycle.
- **Extractor-untouchable (CLAUDE.md)** — this slice IS extractor
  work. Production-discipline: `./smoke_tests` after the slice
  before any 685-fleet rerun.

## 12. Pre-existing WIP (separate from this slice)

- **3 untracked files in `tests/topology/`** from 2026-05-07 (README +
  `test_legacy_amber_semantic_api.cpp` + `test_legacy_amber_semantic_tables.cpp`).
  SDET guardrails on the substrate tables. **Pre-flight attempted
  2026-05-08; tests don't pass against current substrate** — they
  predate `6ec9bff` (Bundle C Slice B ring migration) and `2de56ae`
  (H1 backbone-tag re-bless) which moved fields the tests pin. 1/6
  api failures, 10/14 table failures all look like substrate drift,
  not real bugs. The tracked `topology_semantic_integration_tests`
  IS green on current substrate and is what certifies the substrate
  this slice consumes. Updating the api/table tests is a separate
  forensic pass; out of scope for this slice. Files stay untracked.
- **Stash `agent-integration-WIP`.** Obsolete (the first attempt at
  substrate composition; refined and landed in commits `43c88d4`,
  `ba647cd`, `a5b4732`, `23cf789`). **Do not pop.**
- **Task #10 substrate coverage test.** Pending from the 2026-05-05
  evening handoff. The 30 emitted tables aren't yet programmatically
  sanity-checked beyond build-clean + StringBarrier. Locked decision:
  this slice precedes Task #10.

## 13. Coordination notes

### Memory keyspace mismatch

The two new memories saved 2026-05-08 are in
`-shared-2026Thesis/memory/` (the prior session was launched from
`/shared/2026Thesis/`, one level up from the project root). Sessions
launched per CLAUDE.md convention from `/shared/2026Thesis/nmr-shielding/`
won't auto-load them. The asymmetry rule is captured verbatim in §2 of
this doc, so duplicating is unnecessary for this slice — just launch
the next session from the conventional dir.

### `~/.nmr_tools.toml` machine-local key

Already added 2026-05-08:
```toml
bmrb_atom_nom = "/shared/2026Thesis/nmr-shielding/references/bmrb_data/atom_nom.tbl"
```

### Scale target

Application is **2300 proteins in RefDB** for codex's Stage 1 stats
redo, not the test fixtures. Edge cases at that scale: legacy BMRB
records with pre-Markley naming, multi-chain proteins, modified
residues. Output-side logged-fallback is the safety valve. No
equivalent on the input side — `NamingApplicator` must be exact across
the whole 2300.

## 14. References

- `references/bmrb_data/atom_nom.tbl` (committed in `6c4d13e`).
- `references/ANNOTATED_BIBLIOGRAPHY.md` — bibliography section on the table.
- `src/generated/LegacyAmberSemanticTables.{h,cpp}` — substrate source of truth.
- `src/SemanticEnums.h` — typed enums consumed in the schema.
- `src/MutationDeltaResult.{h,cpp}` — current matchup (target of rewrite at lines 278-330).
- `src/FramePdbEmitter.{h,cpp}` — the structural precedent.
- `src/RuntimeEnvironment.{h,cpp}` — static-accessor precedent for `BmrbAtomNom()`.
- `spec/plan/naming-applicator-architecture-sketch-2026-05-06.md` — input-side architecture (the asymmetric counterpart).
- `spec/plan/session-handoff-20260505-evening.md` — Task #10 detail.
- `~/.claude/projects/-shared-2026Thesis/memory/feedback_naming_input_output_asymmetry.md`
- `~/.claude/projects/-shared-2026Thesis/memory/feedback_bmrb_chemistry_not_raw.md`
