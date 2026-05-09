# Post-topology Doc-Cleanup — 2026-05-09

**Status:** in flight. Consolidates the loose ends identified across
the 2026-05-08 → 09 topology + projection + SDK + bones-retirement
session set. Goal: doc surface honest before forward calculator
construction begins.

## Why this exists

Several large slices landed yesterday and earlier today:

- AMBER substrate + NamingApplicator + Bundle C ring substrate (commits
  `4ba5491` through `6ec9bff`).
- CategoryInfoProjection + MutationDeltaResult typed-identity matchup
  + dia/para components + SDK wrapper (commits `8accdb6`, `e1b5bcc`,
  `d1ad904`).
- Doc bones-retirement, PlanarGeometryResult and
  AIMNet2PolarisabilityResult planned-calculator amendments, Trp-cage +
  synthetic peptide as workbook substrates (commits `08034cf`,
  `8aaa21a`, `aa3c1f1`, `be880f2`).

Each of those slices left small wakes — a WIP-marker that is no longer
WIP, a stale audit doc, a memory file in the wrong keyspace, an untracked
test file that doesn't pass current substrate. Individually trivial;
together a coherence problem if forward work begins on top of them.

This plan captures every item identified across that conversation and
sequences them so nothing gets lost. After all items complete, the
plan retires to `spec/plan/bones/`.

## Items

### Tier 1 — Doc surgeries (small, decisive)

- [ ] **PATTERNS.md §"(WIP) Rule-application architecture for runtime
      atom-name canonicalisation"** (~line 522) — promote from WIP to
      canonical. The trigger conditions ("after the NamingRegistry
      refactor lands and Session E projections exercise the same
      shape") have both landed (`NamingApplicator` in `4ac7d79`;
      `CategoryInfoProjection` in `8accdb6`). Drop the WIP marker;
      write the actual canonical-pattern prose; cross-reference both
      module names.

- [ ] **OBJECT_MODEL.md §"(WIP) Runtime atom-name canonicalisation"**
      (~line 420) — promote. Same trigger conditions met. Replace the
      "section will be filled in once architecture lands" placeholder
      with the canonical description: typed-rule object holds the
      published rule-sets, per-atom transient application maps,
      explicit resolution method documenting cross-rule-set choices.
      Name `NamingApplicator` (input) + `CategoryInfoProjection`
      (output) as the two ends of the symmetric pattern.

- [ ] **PATTERNS.md "The naming boundary"** (~line 396) — currently
      describes the *input* side only (NamingRegistry as the single
      translation point). Add a paragraph naming the output side —
      `CategoryInfoProjection` is the symmetric counterpart at the NPY
      emission boundary; output names are tolerant where input is
      strict; logged-fallback acceptable; the asymmetry is
      architectural per `feedback_naming_input_output_asymmetry`.

- [ ] **EVIL_STRING_AUDIT_2026-04-28.md** retire to
      `spec/plan/bones/`. Per the 2026-05-09 Explore-agent refresh,
      the surface has resolved: 13 critical correctness sites
      migrated; 8 cosmetic + 1 deferred chi-angle slice remain. Add a
      final-state footer summarising the resolution before moving the
      doc.

- [ ] **KNOWN_BUGS.md** — audit pass. The H1-backbone-tag bug noted
      there was fixed in commit `85de965` plus re-bless in `2de56ae`.
      Re-read each entry against current `master` (`be880f2`); remove
      or mark-resolved items the topology + projection slices closed;
      confirm the chi-angle resolver deferral is documented as the
      remaining string site. Was last touched pre-2026-05-04 — likely
      stale on multiple entries.

### Tier 2 — Decisions + work (medium)

- [ ] **3 untracked `tests/topology/*` files** — README +
      `test_legacy_amber_semantic_api.cpp` +
      `test_legacy_amber_semantic_tables.cpp`, all from 2026-05-07.
      Decision: **retire to `spec/plan/bones/`** (or
      `tests/bones/topology/`) with an archaeology note. The api and
      tables tests' assertions captured a useful substrate-table
      shape, but they predate Bundle C Slice B (commit `6ec9bff`,
      ring migration to `LegacyAmberTopology`) and 10/14 + 1/6
      assertions fail against current substrate. Re-writing the
      assertions for current substrate is larger than this cleanup
      pass; the tracked `topology_semantic_integration_tests`
      provides live coverage. Bones'd files are preserved for future
      substrate-aware re-write if anyone wants the audit shape back.

- [ ] **Memory keyspace** — copy
      `feedback_naming_input_output_asymmetry.md` and
      `feedback_bmrb_chemistry_not_raw.md` from
      `~/.claude/projects/-shared-2026Thesis/memory/` (parent-dir
      keyspace, where the 2026-05-08 crashed session was running) to
      the canonical `~/.claude/projects/-shared-2026Thesis-nmr-shielding/memory/`
      keyspace per CLAUDE.md "Working-directory convention." Update
      that keyspace's `MEMORY.md` index. Both memories codify durable
      architectural rules that should auto-load from sessions started
      conventionally.

- [x] **NMR_EXTRACT_DESIDERATA_2026-04-22.md Done-vs-Pending pass**
      — Status header landed 2026-05-09. Per the Explore-agent audit,
      5 items LANDED, 3 PARTIAL, ~25 PENDING. The doc does NOT retire
      to bones/ in this pass — for ~25 PENDING items it remains the
      only canonical design-intent home. Doc retirement requires
      first migrating those items to canonical pending homes
      (PLANNED_CALCULATORS amendments, comprehensive-inventory
      amendments, or new I/O / diagnostics docs). That migration is
      Tier 2.5 (separate ~2–3 hour slice).

- [ ] **Comprehensive inventory Section 10 IMPROVEMENT entries** —
      three items: ring-normal stability fix
      (BiotSavartResult + HaighMallionResult), single-loop default
      audit (BiotSavartResult), H-bond geometry angle θ over distance
      d. Verify whether any landed during the 2026-05-08 → 09 work
      and update status accordingly.

### Tier 2.5 — DESIDERATA item migration (separate slice, deferred)

Surfaced by the 2026-05-09 Explore-agent audit. ~25 items in
`spec/NMR_EXTRACT_DESIDERATA_2026-04-22.md` Sections A–E have design
intent (physics references, dependencies, validation-bench mapping)
recorded only in that document. Migrating them to canonical pending
homes is required before DESIDERATA can retire to bones/. Estimated
~2–3 hours; not in the original cleanup pass.

Items needing migration, by destination:

- **PLANNED_CALCULATORS_2026-04-22.md amendments** (new
  ConformationResults: A.1 CSAPrincipalAxis, A.2 AmideTensorGeometry,
  A.3 BulkSusceptibilityAccumulator, A.4 NICSProbe,
  A.6 ParamagneticRelaxationEnhancement, A.8 MemoryKernelExtraction,
  A.9 ErgodicityMetric, A.10 CCRRate, A.11 BenchmarkBackCalculation).
  ~9 amendments.
- **comprehensive-calculator-inventory-2026-04-30.md Section 10
  amendments** (kernel-pass improvements: B.2 smooth cutoffs,
  B.3 multipole Coulomb, B.4 distributed ring current,
  B.5 H-bond cooperativity, B.6 MOPAC shielding_contribution,
  B.7 C8 dispersion). 6 entries; some may already overlap with
  Section 10 entries — verify and consolidate.
- **New `spec/I_O_AND_SCHEMA_2026-XX-YY.md`** consolidating Section C
  items (C.1–C.8: ORCA bond orders, ExperimentalReferenceLoader,
  lanthanide inputs, magnetizability, irrep metadata, units/signs,
  bench H5 slices, TrajectoryResult schema formalisation).
- **New `spec/DIAGNOSTICS_AND_WORKFLOWS_2026-XX-YY.md`** consolidating
  Section D and applicable Section E items (D.2 T2-residual map,
  D.3 T2-correlation map, E.2 event menu, E.3 probe-point API,
  E.6 per-bench pipeline mode). E.5 calibration-raw enforcement
  goes in PATTERNS.md as a discipline rule.

After migration completes, DESIDERATA retires to `spec/plan/bones/`
with a final-state footer like the EVIL_STRING audit.

- [ ] Tier 2.5 — DESIDERATA item migration (deferred to its own slice)

### Tier 3 — Real engineering work (delicate)

- [x] **Chi-angle resolver migration** in
      `CacheResidueBackboneIndices_Typed` (`Protein.cpp:728-751`)
      complete 2026-05-09. Each chi-atom name from
      `AminoAcidType::chi_angles` is now converted to an
      `AtomMechanicalIdentity` via
      `topology_generated::ComputeAtomMechanicalIdentity` and looked
      up via `LegacyAmberTopology::ResidueAtomsWithIdentity`. Strings
      stay at the AAType-table boundary; the lookup is typed.
      Bit-equivalent for canonical AMBER input on this fixture. The
      first-pass chi-angle string match in
      `CacheResidueBackboneIndices` (lines 624-634) is left in place
      as the loader-boundary fallback — that's sanctioned per
      PATTERNS.md "Naming boundary." Tests:
      LegacyAmberSemanticIntegration 14/14, CategoryInfoProjection
      10/10, MutationDelta 14/14, structure-tests chi/DSSP filter
      13/13.

- [ ] (DEFERRED, codex-conditional) **`pseudoatom_super_group` enum**
      on the structured NPY — currently emit only `in_super_group`
      bool; codex derives QG/QD/QH/QR from
      `(residue_type, locant, in_super_group)`. Promote to typed
      column only if codex asks during the Stage 1 stats redo. Don't
      add speculatively.

- [ ] (DEFERRED, codex-conditional) **`match_kind` enum** on
      `MutationDeltaResult` per-atom NPY — currently aggregate
      counts in `OperationLog::Info` only. Promote to per-atom
      column if codex needs row-level rejection categorization.
      Don't add speculatively.

### Tier 4 — Pre-flight for next forward slice

- [ ] **Trp-cage + synthetic peptide PDBs** at
      `tests/data/illustrative_peptides/` per the
      `PLANNED_CALCULATORS_2026-04-22.md` Amendment 2026-05-08(c)
      decision. Build via `tleap`:
      - `1l2y_trpcage.pdb` from PDB ID 1L2Y first model (or load via
        cifpp + extract a single conformation)
      - `synthetic_22mer.pdb` from sequence `ACDEFGHIKLMNPQRSTVWYCV`
        with `tleap` extended-conformation builder
      - Disulfide bond between the two Cys residues in the synthetic
        (positions 2 and 21 in the sequence)
      Commit alongside a `BUILD.md` listing the tleap commands so
      regeneration is deterministic. ~5 min tleap.

- [ ] **AIMNet2 `.jpt` requires_grad pre-flight check** per
      `PLANNED_CALCULATORS_2026-04-22.md` Amendment 2026-05-08(b).
      The 10-line Python script: load
      `data/models/aimnet2_wb97m_0.jpt`, run forward with
      `requires_grad=True` on coordinates, backward on `sum(charges)`,
      assert `coords.grad` propagates. Confirms (or refutes) whether
      `AIMNet2PolarisabilityResult` can be built without re-exporting
      the model from PyTorch. Gates the calculator slice. ~10 min.

## Acceptance criteria

- All Tier 1 items committed (one commit per logical group is fine,
  or one combined doc-cleanup commit; each commit message names the
  specific changes).
- Tier 2 decisions executed per the recommendations above; the
  judgment calls (retire-not-fix for `tests/topology/*`, copy memory
  to canonical keyspace) documented in commit messages.
- Tier 3 chi-angle migration lands as its own commit with targeted
  test verification before/after. The deferred codex-conditional
  items stay deferred — do not add speculatively.
- Tier 4 PDBs committed with deterministic-regeneration BUILD.md;
  pre-flight Python check result captured (in this doc's completion
  footer) before retirement.
- This plan doc itself retires to `spec/plan/bones/` after the final
  item lands, with a footer summarising what was done.

## Order of execution

1. Tier 1 (~60 min, single coherent commit) — pure doc surgery.
2. Tier 2 (~60–90 min, may be 2–4 commits) — decisions + work.
3. Tier 4 PDBs (~15 min, one commit with BUILD.md).
4. Tier 4 AIMNet2 pre-flight (~15 min, results documented in the
   completion footer of this plan).
5. Tier 3 chi-angle migration (~60 min, its own commit with test
   verification).
6. Plan retirement (~5 min) — move this doc to
   `spec/plan/bones/`, then commit.

Total estimated cost: 3–4 hours across the work.

## On completion

The plan retires to `spec/plan/bones/` with the AIMNet2 pre-flight
result and any decisions made during execution preserved as a
footer.
