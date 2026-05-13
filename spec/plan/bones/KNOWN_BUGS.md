# Known bugs and architectural debt

## 2026-05-09 update — substantially all resolved

The 2026-04-29 AMBER charge slice + 2026-05-06 NamingApplicator +
2026-05-07 Bundle C ring substrate + 2026-05-08 CategoryInfoProjection
+ MutationDeltaResult typed-identity matchup landings, taken together,
**resolve every bug originally documented here**. Status table
matched against current `master`:

```text
Bug 1 (silent ff14SB charge fallback)          RESOLVED. ParamFileChargeSource
                                               loud-fails on missing terminal
                                               templates and missing rows.
                                               Six tests in
                                               test_foundation_results.cpp.
                                               Fixed in the AMBER charge slice
                                               (commits 4ba5491-91e4124).

Bug 2 (string-keyed dispatch ~50 sites)        RESOLVED. Per the 2026-05-09
                                               EVIL_STRING audit refresh
                                               (now retired to
                                               spec/plan/bones/), all 13
                                               critical-correctness sites
                                               migrated to typed-substrate
                                               access. 8 cosmetic sites
                                               (PDB formatting, error
                                               messages, display labels)
                                               correctly retain string use
                                               at human-facing surfaces.
                                               1 chi-angle resolver site
                                               deferred — tracked in
                                               spec/plan/post-topology-doc-
                                               cleanup-2026-05-09.md Tier 3.

Bug 3 (deferred CHARMM↔IUPAC rules)           RESOLVED. The
                                               commented-out rules
                                               (NamingRegistry.cpp:172-307
                                               originally) are now active in
                                               the NamingApplicator (commit
                                               4ac7d79 and follow-ups).
                                               GLY pre-Markley, PRO/LYS/ARG
                                               δ-methylene, ILE δ-methyl,
                                               ALA β-methyl wildcard fix
                                               all in place per
                                               NamingRegistry.cpp:250-289.

Bug 4 (pdb_atom_name forgotten-source string)  RESOLVED. All 4 loaders
                                               canonicalise through
                                               GlobalNamingApplicator() at
                                               the boundary. The
                                               LegacyAmberTopology +
                                               ForceFieldChargeTable
                                               landing typed the implicit-
                                               AMBER assumption rather than
                                               enforcing it. Atom names
                                               surface as typed substrate
                                               fields downstream
                                               (SemanticAt(ai).backbone_role,
                                               .locant, etc.); the
                                               pdb_atom_name string remains
                                               on Atom for output
                                               projection only.

Design choice 1 (KD-tree position matching)    RESOLVED. MutationDeltaResult
                                               typed-identity matchup
                                               landed in commit e1b5bcc
                                               (2026-05-08): binds by
                                               (residue_index,
                                               AtomMechanicalIdentity);
                                               consume-in-order for
                                               equivalent-H sets;
                                               spatial-NN diagnostic on
                                               rejections. KD-tree +
                                               nanoflann + MATCH_TOLERANCE
                                               removed.

Regression 1 part 1 (delta channel collapse)   RESOLVED. Six new dia/para
                                               NPYs landed in commit
                                               e1b5bcc: wt_/mut_/delta_ ×
                                               diamagnetic/paramagnetic.
                                               σ_total = σ_dia + σ_para
                                               within ORCA output precision.

Regression 1 part 2 (ORCA load silent-Ok)      PARTIAL. The
                                               OperationRunner per-conformation
                                               failure path propagates correctly
                                               (RunMutantComparison checks
                                               WT/ALA wt_result.Ok at
                                               OperationRunner.cpp:251 and
                                               260). However, ORCA
                                               shielding-load failure
                                               specifically still logs Error
                                               and continues silently
                                               (OperationRunner.cpp:218-223).
                                               This is a small focused fix
                                               separate from the
                                               typed-substrate work — change
                                               line 222 to set out.error and
                                               return early. Tracked as a
                                               near-term fix but not in the
                                               post-topology cleanup pass
                                               since it is calculator-side
                                               correctness, not doc surgery.

Regression 1 part 3 (collision-fatal in        RESOLVED-BY-DESIGN. The
                       cross-protein lookup)   typed-identity consume-in-order
                                               matcher binds each mut atom
                                               at most once by construction.
                                               No collision check is needed
                                               because collisions cannot
                                               occur — equivalent-H sets
                                               within a residue are bound
                                               position-by-position from a
                                               consumed pool.
```

The pointer to the active plan section at the bottom of this file
(`spec/plan/openai-5.5-strong-architecture-layout.md`) remains the
architectural reference. Operational reference is memory entry
`project_proteintopology_architecture`.

The rest of this document is preserved as the original honest list of
bugs as documented post-revert. The descriptions remain diagnostic and
historically accurate; the "deferred indefinitely" framing is overtaken
by the resolution table above. **One residual remains — Regression 1
part 2 (ORCA silent-Ok at `OperationRunner.cpp:218-223`).**

---

## Original entry (pre-2026-04-29 framing)

**Scope.** Issues known to exist in the working tree at commit `130c8da`
(post-revert state, 2026-04-27). Documented honestly rather than fixed —
see "Why this document exists" below.

**Why this document exists.** A 2026-04-26 attempt to add a typed IUPAC
topology layer to `Atom` introduced architectural problems that would
have required a 10–20 session rebuild to resolve cleanly. Given finite
two-person-team capacity (the user + Claude), the rebuild's opportunity
cost (weeks of blocked thesis work, equipment costs accumulating) was
larger than its benefit. The decision was to revert the IUPAC layer and
document the bugs it surfaced — so the working tool remains productive
while the bugs are known and can be addressed when resources permit.

The active plan that emerged from the post-revert design conversation
is the **`LegacyAmberTopology` + calculator-contract** architecture,
captured authoritatively in
`spec/plan/openai-5.5-strong-architecture-layout.md` and operationally in
memory entry `project_proteintopology_architecture` (loaded
automatically at session start). The IUPACAnnotation framing was a
stage in the design; the typed semantic fields it would have provided
are absorbed into `LegacyAmberTopology`. See "Pointer to the active
plan" at the bottom of this file. The full design evolution (six
framings between 2026-04-26 and 2026-04-28) is recorded in the
memory entry's Historical section.

---

## Bug 1 — Silent charge corruption on `--pdb` / `--protonated-pdb` / `--orca` paths

**Severity.** Physics-corrupting on affected paths. Trajectory path
unaffected.

**Where.** `src/ChargeSource.cpp:48` (and the fallback at line 56).

**What happens.** `ParamFileChargeSource::LoadCharges` builds an ff14SB
parameter-file lookup key as:

```
key = ff_resname + " " + identity.pdb_atom_name
```

The ff14SB parameter file is keyed by AMBER atom names (CA, CB, HB1,
HB2, etc.). When `pdb_atom_name` is not exactly the AMBER convention
(e.g., a CHARMM-named PDB supplied to `--protonated-pdb`, or a non-AMBER
ORCA path), the primary key misses. The fallback at line 56 tries the
canonical residue name with the same atom name; that also misses.
Control flow falls through to `ChargeSource.cpp:62-72`, which silently
assigns element-default charges (`charge = 0.0`, radius from a
hard-coded element table) with no warning, no diagnostic, no log.

**Affected load paths:**
- `--pdb` (reduce-protonated PDB → AMBER): mostly OK because reduce
  emits AMBER-convention names; can break if the input PDB had
  non-standard naming that reduce passed through.
- `--protonated-pdb` (user-attested already-protonated): broken if the
  user's input PDB used CHARMM, GROMACS-source, or other non-AMBER
  naming.
- `--orca` (tleap-built prmtop): broken if the prmtop's ATOM_NAME entries
  diverge from ff14SB-table conventions for any atom.
- `--mutant` (uses ORCA twice): same as `--orca`.

**Unaffected:**
- `--trajectory` (full-system XTC + TPR): the loader translates CHARMM →
  Standard via `NamingRegistry::TranslateAtomName` at
  `FullSystemReader.cpp:705`, so the lookup keys are correct.

**Operational guidance.** For production extractions where charge
correctness matters (calibration runs, Coulomb-derived shielding,
EFG calculations), prefer the `--trajectory` path. Or audit the
specific input proteins: the bug fires when `(ff_resname,
pdb_atom_name)` is not in `data/ff14sb_params.dat`.

**Detection.** No diagnostic surface in the current code. To detect, post
an extraction with `aimnet2_charges.npy` (or any charge-emitting NPY)
and check for atoms with `partial_charge == 0.0` whose element wouldn't
naturally produce zero.

**Architectural debt context.** This bug was identified during the
2026-04-27 audit work and is pre-existing — it existed before the
IUPAC topology landing and survives the revert.

---

## Bug 2 — String-keyed dispatch on atom names throughout `src/`

**Severity.** Architectural debt, not a runtime correctness bug in the
current code. Becomes a correctness bug in any future code that violates
the implicit "atom names are AMBER-convention" assumption.

**Where.** Approximately 50 sites across `src/` per the audit's Phase 1
inventory. Representative sites:

- `src/Protein.cpp` — `CacheResidueBackboneIndices` matches `pdb_atom_name`
  against literals "N", "CA", "C", "O", "H", "HN", "HA", "HA2", "CB".
- `src/Protein.cpp` — `DetectAromaticRings` builds a string-keyed map
  per residue and looks up ring atom names from `AminoAcidType.rings[k]`
  (currently `vector<const char*>`).
- `src/Protein.cpp` — HIS tautomer detection by string match on "HD1"
  / "HE2".
- `src/ProtonationDetectionResult.cpp` — variant detection by string
  matches on "HD1", "HE2", "HD2", "HG", "HZ1", "HZ2", "HZ3".
- `src/CovalentTopology.cpp` — bond classification keyed in places by
  atom name strings.

**Why this is debt rather than a bug — historically.** The system has
worked because every load path either (a) translates to AMBER
convention via `NamingRegistry::TranslateAtomName` (full-system
trajectory path) or (b) happens to be supplied AMBER-convention input
(reduce / tleap on the non-trajectory paths). The string-keyed
dispatch is a load-time boundary operation and the bounded universe
of inputs has historically been AMBER-convention.

**Why it is not actually debt — the structural critique.** See Bug 4.
The `pdb_atom_name` field is a forgotten-source string with
inconsistent translation discipline across loaders, read by ~50
consumer sites that all assume AMBER convention. The "happens to be
AMBER-convention input" assumption is not enforced at any boundary;
it is documentation. Bug 1 (silent ChargeSource corruption) is the
acute symptom. The full-coverage fix is scheduled (see Bug 4), and
runs **before** the IUPACAnnotation work, not after.

**Architectural debt context.** Identified during the 2026-04-27 audit
work; full inventory preserved with the audit material in the
emergency archive (do not extract).

---

## Bug 4 — `pdb_atom_name` is a forgotten-source string

**Severity.** Architectural correctness bug. The acute symptom is Bug
1; the architectural shape encompasses Bug 2's full ~50-site dispatch
surface plus the loader-side inconsistency.

**Where.** `Atom::pdb_atom_name` is a `std::string` field on every
Atom. Four loaders write to it. ~50 sites across `src/` read it.

**The structural problem.** The field's name implies "the atom's
name from the PDB." The reality is "a string from one of four
sources, with inconsistent translation discipline, read by consumers
that assume AMBER convention regardless of source."

| Loader | Translation discipline |
|---|---|
| `FullSystemReader` | NamingRegistry::TranslateAtomName → AMBER (line 705) |
| `GromacsEnsembleLoader` | NamingRegistry::TranslateAtomName → AMBER (lines 337, 517) |
| `OrcaRunLoader` | Residues canonicalised via NamingRegistry; **atom names assigned raw** (line 204) |
| `PdbFileReader` | Raw cif++ label, **no translation** (line 131) |

Downstream consumers (`ChargeSource` ff14SB lookup,
`Protein::CacheResidueBackboneIndices`, `Protein::DetectAromaticRings`,
`ProtonationDetectionResult` variant dispatch, HIS tautomer
detection, `CovalentTopology` bond classification, ~50 sites total)
read `atom.pdb_atom_name` and assume AMBER convention. When the
loader didn't translate and the input wasn't AMBER-convention to
begin with, every consumer becomes a silent-corruption site.

**Why this is the bug shape, not "debt."** A field whose semantics
depend on which loader produced it, with no compile-time or runtime
type that records which loader did, is broken by construction. The
"singular-atom collapse" framing the user has articulated repeatedly:
*atoms are not named — there are strings from sources, and that's
typed work.* The current field pretends otherwise.

**Why this is dissolved by `LegacyAmberTopology`, not patched in
isolation.** A previous framing proposed fixing the boundary discipline
as a standalone session before IUPACAnnotation. Subsequent design
review (2026-04-28) identified that fix as cementing AMBER-as-implicit-
truth: every validator-based fix declares "the protein must be AMBER,"
which makes AMBER the system's canonical force field by code rather
than by choice. Adding any new layer on top would inherit that.
The dissolution comes from typing the contract itself: `LegacyAmberTopology`
is the *named* AMBER calculator contract, not the implicit one;
loaders produce instances of it under explicit translation; calculators
declare they require it; future contracts (CHARMM peer, IUPAC peer)
plug in as peers without privileging the legacy.

**Resolution shape.** Captured in
`spec/plan/openai-5.5-strong-architecture-layout.md` and the
`project_proteintopology_architecture` memory entry. Summary:

1. `Atom::pdb_atom_name` becomes a `LegacyAmberTopology`-owned typed
   value, not a free string. Loaders translate input strings to AMBER
   convention at the boundary (PdbFileReader and OrcaRunLoader gain
   `NamingRegistry::TranslateAtomName` calls, with `AminoAcidType`
   validation as the post-check); CHARMM input via GROMACS uses the
   existing CHARMM↔Standard rules. Non-AMBER atoms fail loudly at the
   loader, not silently at consumers.
2. The 28 dispatch sites move off raw `pdb_atom_name` reads onto
   `LegacyAmberTopology`-typed accessors (typed slots, typed enrichment
   flags, typed protonation variant). The typed accessors are the
   contract; the AMBER strings are projection of typed values back to
   convention strings at output time.
3. `ChargeSource` string-dispatch (the silent-0.0 surface in Bug 1)
   moves to loader-time inside `LegacyAmberTopology` resolution.
   Missing keys produce loud errors; `ParamFileChargeSource`'s
   silent-fallback path is deleted.
4. Test coverage of all four loaders × representative dispatch sites ×
   non-AMBER inputs. NPY parity at every commit boundary throughout
   the calculator sweep.

**Architectural debt context.** Identified 2026-04-28 as the bug shape
underlying Bugs 1 and 2. Dissolved by the `LegacyAmberTopology`
landing rather than fixed in isolation. Pre-spec calculator × topology
matrix work runs first to bound the sweep into single-session
migrations.

---

## Bug 3 — Deferred `NamingRegistry` rules for CHARMM↔IUPAC gaps

**Severity.** Affects atoms produced by the `--trajectory` path on
specific positions in specific residues. Not a runtime crash; produces
wrong labels in the H5 output for those positions, and (less obviously)
mis-assigns the atoms during boundary translation if the wildcard rule
fires incorrectly.

**Where.** `src/NamingRegistry.cpp:172-307` (the commented-out rule
block).

**What's missing.** Eight categories of CHARMM↔IUPAC translation rules
were left commented out pending fleet-wide vetting. They cover:

- ILE γ-carbon, δ-methyl, γ1-methylene atom name translations.
- δ-methylene rules for ARG, LYS, PRO.
- ε-methylene rule for LYS.
- α-methylene rule for GLY.
- ALA β-methyl identity blockers (the wildcard β-methylene rule
  over-fires on ALA, producing duplicate `HB3` labels).

**Operational guidance.** For 10-protein calibration work, two typed
Python data constants in `h5-reader/notes/nmr_forensics/pack_experimental_shifts.py`
handle the affected cases in-place (per the comments in
`NamingRegistry.cpp`). Fleet-wide vetting is needed before activating
the rules in the C++ registry.

**Architectural debt context.** Documented in
`spec/ChangesRequiredBeforeProductionH5Run.md` on master.

---

## Design choice 1 — `MutationDeltaResult` uses KD-tree position matching

**Severity.** Not a bug. A known design choice with a planned resolution.

**Where.** `src/MutationDeltaResult.cpp` (pre-IUPAC implementation,
which is what the current tree has).

**What it does.** Cross-protein WT-to-mutant atom matching is done by
KD-tree on positions with a 0.5Å tolerance. This is the matcher that
produced the existing 723-protein calibration corpus.

**Why this is a design choice rather than a bug.** It works. The
0.5Å tolerance is tight enough to avoid false matches in practice; the
calibration results are stable.

**Resolution path (planned, not yet built).** Under the
`LegacyAmberTopology` plan (`spec/plan/openai-5.5-strong-architecture-layout.md`
+ memory entry `project_proteintopology_architecture`),
`MutationDeltaResult` declares `Contract = CalculatorContract<LegacyAmberTopology, ...>`
and migrates from KD-tree position matching to typed
`(residue_index, atom-slot-via-AminoAcidType, element)` matching using
the LegacyAmberTopology accessors. No more geometric tolerance;
cross-protein matching becomes typed-by-chemistry. Targeted scope
(~50-line change in MutationDeltaResult.cpp once LegacyAmberTopology
is in place).

---

## Regression 1 — MutationDelta channel collapse + run-boundary silent-success

**Severity.** Two regressions introduced by the 2026-04-27 revert.
Bundled into one entry because they ride together — both will land
alongside the IUPACAnnotation typed-matching migration above.

**What we briefly had (commits `7ce1d86` + `702f794`, now reverted):**

1. `MutationDeltaResult` emitted nine tensors per matched atom — WT,
   mutant, and delta × {total, diamagnetic, paramagnetic} — instead
   of the single `delta_shielding_total` channel. The reasoning
   recorded in the commit: paramagnetic channel is where heavy-atom
   and aromatic-ring effects dominate; residual analysis on the total
   alone cannot distinguish dia vs para. Constitution-level
   anti-simplification.
2. `OperationRunner` hard-failed on `--orca` load failure (was
   silently returning Ok with the result absent — a 723-protein batch
   with a typo'd ORCA path produced WT-only outputs that look like
   success). And on `RunMutantComparison` MutationDelta failure (was
   logging Error but not propagating).
3. Cross-protein matching enforced collision-is-fatal: two atoms in
   the same residue with the same IUPAC name surfaced loud instead
   of being silently de-duplicated.

**Why deferred rather than cherry-picked.** User direction 2026-04-28:
"mutation delta and anything of that nature we want. It will have to
wait for the annotation though, the way we associate atoms without
residue and info in the comparison is not worthy." The current
post-revert KD-tree-on-positions matcher is the unworthy one;
re-landing per-channel decomposition + hard-fail + collision-fatal on
top of it would marry good discipline to a matching scheme we plan to
replace. Land them all together when MutationDeltaResult migrates to
typed matching under `LegacyAmberTopology`.

**Current behaviour (post-revert):**
- `MutationDeltaResult` emits delta_shielding_total only; dia and para
  channels are computed by ORCA and stored on ConformationAtom but
  not differenced.
- `OperationRunner` with bad `--orca` path silently returns Ok with
  no shielding result; a `RunMutantComparison` with a missing delta
  logs Error and returns Ok.
- Cross-protein matcher: KD-tree at 0.5Å tolerance, no collision check.

**Resolution.** Bundled with the `LegacyAmberTopology` MutationDelta
migration. When that work opens, add to its scope:
- 9-tensor channel decomposition (WT/mut/delta × total/dia/para)
- OperationRunner hard-fail on ORCA load + MutationDelta failure
- Collision-fatal in the typed cross-protein lookup map

---

## Pointer to the active plan

The active architecture is **`spec/plan/openai-5.5-strong-architecture-layout.md`**
(authoritative). The operational reference is memory entry
`project_proteintopology_architecture` in the user's persistent memory
store, loaded automatically at session start.

In one paragraph: `LegacyAmberTopology` is the explicit typed contract
for current calculators — residue symbolic slots (N/CA/C/O/H/HA/CB,
chi tuples, protonation variant index), ring topology with all
type-specific parameters (Intensity, JBLobeOffset, NitrogenCount,
etc. — ring typing changes calculator math), bond topology
(`CovalentTopology` becomes the bond-graph component inside),
legacy enrichment projection (AtomRole, Hybridisation, is_backbone,
is_amide_H, etc.), and the legacy AMBER atom-naming as typed values
on the contract. `Atom::pdb_atom_name` becomes a contract-owned
typed value. Calculators declare
`using Contract = CalculatorContract<LegacyAmberTopology, ChargeSetT>`;
template helpers (`RequiredTopology<CalculatorT>(conf)`,
`RequiredCharges<CalculatorT>(conf)`) route compatibility checks at
instantiation. `ForceFieldChargeSet` is a first-class object on
Protein with typed source kind (AmberPrmtop, GromacsTpr,
Ff14SBParamFile, Preloaded). Calculated charges remain result-owned.
Naming conventions (IUPAC, BMRB) are output projections of
`LegacyAmberTopology`, not separate contracts. Future richer topology
contracts plug in as peers — `LegacyAmberTopology` is "just another
calculator topology, not a special one." Matrix pre-spec runs first
(calculator × topology-field matrix); round-robin design with two
Opus + one OpenAI 5.5, agreement per stage.

The earlier framings — IUPAC topology layer on Atom (2026-04-26,
reverted), TopologyProtein ABC (2026-04-27, scaled back as version-2-
shaped), IUPACAnnotation as targeted minimum (2026-04-28 morning),
evil-string-fix-first (2026-04-28 mid-day), typed FF + typed
convention (2026-04-28 afternoon) — are recorded as design history in
the memory entry's Historical section. The 2026-04-26 IUPAC working
state is archived to
`/shared/2026Thesis/iupac-fix-attempt-archive-2026-04-27.tar.gz`
(~275MB, outside the repo, **do not extract**). The memory entry is
self-contained; the archive should not need extraction.
