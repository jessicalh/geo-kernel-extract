# Session handoff — 2026-04-30 evening

**Purpose:** Pickup brief for the morning continuation. Read this first; it points at everything else.

**Written at end of long session 2026-04-30 (poster timeline 8 weeks, results 14 weeks). User's intent for the morning: review the deferred OpenAI critiques, refine the remaining topology-capture question (CYS / HIS / ASH / GLH / LYN real-chemistry encoding), then start landing NC1-NC7 implementation.**

---

## State of the working tree

```
HEAD (origin/master)        f8729e3  KNOWN_BUGS: align with IUPACAnnotation plan
                            (from 2026-04-26)

UNCOMMITTED — AMBER charge slice (six steps; landed 2026-04-29 evening)
  ~12 new src files + ~20 modified, ~2000 lines C++
  62/62 acceptance tests passing
  Per-step commit boundaries intentional per the central plan
  STATUS: commit-ready, awaiting deliberate per-step commits
  Spec:   spec/plan/amber-implementation-plan-2026-04-29.md

UNCOMMITTED — PHASE 0 partial (started 2026-04-30 morning)
  tests/test_amber_trajectory.cpp   (162 lines, untracked)
  src/FullSystemReader.cpp          (+14 lines, modified)
  Smoke test passed (libgromacs reads AMBER TPR cleanly)
  BuildProtein hit "unknown residue HISH" → triggered the wander
  STATUS: smoke + heavy fixture path in tree; loader-side resolution
            is the NC1-NC7 implementation work

UNCOMMITTED — Naming canonical vocabulary spec
  spec/plan/naming-canonical-vocabulary-2026-04-30.md   (33 KB / 774 lines)
  Resolved HISH architecturally (read upstream authorities; do not
    re-decide naming on the consume side)
  Implementation steps NC1-NC7 defined
  STATUS: spec ready; no implementation yet

UNCOMMITTED — All today's planning docs
  spec/plan/running_plan_notes.md
  spec/plan/thesis-discussion-part1-2026-04-30.md
  spec/plan/md-rerun-685-discussion-priors-2026-04-30.md
  spec/plan/md-session-handoff-2026-04-30.md
  spec/plan/comprehensive-calculator-inventory-2026-04-30.md
  (plus this file)

UNCOMMITTED — Spec doc revisions during the day's planning passes
  CLAUDE.md (fleet state corrected)
  spec/TRAJECTORY_RESULT_PLAN_2026-04-24.md (485-cancelled)
  spec/IDENTITY_AND_DYNAMICS_ROLLUP_2026-04-22.md (485-cancelled)
  spec/INDEX.md (modified)
  spec/TEST_FRAMEWORK.md (modified)
  spec/DEPENDENCIES.md (modified)
  KNOWN_BUGS.md (modified)
  OBJECT_MODEL.md (modified)
```

The working tree is in a "lots of completed work, none committed" state. This is fine — the slice was deliberately held at commit-boundary discipline; the rest accumulated during the day's planning. The morning's first decision is whether to land commits before further work or after the topology refinement.

---

## How HISH was resolved (the wander's punchline)

Tried initially: add `HISH` / `HISE` / `HISD` aliases to `NamingRegistry`. User stopped this immediately.

Realised: the prep tree already contains the upstream authorities that record the canonical chemistry. Specifically:

```
prep_run_*/topol.top                     Authority A
                                         The "; residue N <name> rtp <rtp> q <q>"
                                         comment line records the canonical
                                         AMBER name (rtp field = HID / HIE /
                                         HIP / CYX / ASH / GLH / LYN / CYM)
                                         even when the .name field carries
                                         the FF-port label (HISH, HISE, HISD,
                                         CYS-with-HG-missing, etc.).

prep_run_*/sources/                      Authority B
  structure_extracted_*.pdb              The deposit-canonical PDB (or
                                         pdb_companion_*.pdb if structure_
                                         extract didn't run). Atom names
                                         per CCD/IUPAC for the per-atom
                                         canonicalisation step.
```

**The consume side reads, walks, audits, applies. It never re-decides chemistry.** Disulfide pairing was decided by GROMACS specbond.dat at prep time; the rtp record is GROMACS's typed self-report of that decision. Histidine protonation was decided at prep time (interactive `-his` choice or geometric classifier); the rtp record holds it. The trajectory loader reads these records and applies them.

This shape is the canonical-vocabulary spec. Implementation = NC1-NC7.

---

## Real-chemistry capture: use existing typed substrate fields

User asked end-of-day what we'd decided about capturing CYS / HIS / ASH / GLH / LYN real chemistry encoded in the GROMACS intermediate format. Recovered from the central plan (`amber-implementation-plan-2026-04-29.md`):

**The typed substrate already has the fields. We use the existing ones; we do NOT add new ones.**

```
Existing typed encoding the trajectory loader populates:

  Residue.protonation_variant_index           HIS → HID/HIE/HIP via the
                                              AminoAcidType variants table.
                                              CYS → CYS or CYX. ASP → ASP
                                              or ASH. GLU → GLU or GLH.
                                              LYS → LYS or LYN. Etc.

  Residue.terminal_state                      N-/C-terminal variant
                                              (under AmberPreparationPolicy
                                              for the prepared-charge path;
                                              read directly for trajectory
                                              loads).

  CovalentTopology bond list,                 The SG-SG inter-residue bond
    filtered by BondCategory::Disulfide       IS the disulfide record. No
                                              separate disulfide_partners
                                              field on Residue. Calculators
                                              that need "is this SG part of
                                              a disulfide?" filter the bond
                                              list.

  ConformationAtom.partial_charge             Already varies per variant
                                              via ChargeAssignmentResult +
                                              AminoAcidType variants table.
                                              SG charge is -0.08 e for CYX,
                                              -0.31 e for CYS, by lookup
                                              not by inference.

  Atom-set deltas                             IMPLICIT in the typed Protein.
                                              The typed Protein only contains
                                              atoms that ARE present in the
                                              source. CYX has no HG
                                              ConformationAtom because
                                              pdb2gmx didn't write one. HIP
                                              has both HD1 and HE2 because
                                              pdb2gmx wrote both.
```

**The translating we do** is populating these existing typed fields from the upstream-authority files:

```
topol.top rtp comment line          → Residue.protonation_variant_index
                                       (rtp HIP → variant index for HIP;
                                        rtp CYX → variant index for CYX;
                                        rtp ASH → variant index for ASH; etc.)

topol.top [bonds] section + rtp=CYX → BondCategory::Disulfide entry in
                                       CovalentTopology
                                       (each SG-SG bond between two CYX
                                        residues becomes one Disulfide bond
                                        in the typed list)

ChargeAssignmentResult              → ConformationAtom.partial_charge
   (existing pipeline)                from AminoAcidType variants lookup
                                       given residue.type +
                                       residue.protonation_variant_index

Atom-set                            → ConformationAtom records
   (existing — what's present is      Whatever atoms the trajectory
    what's there)                      contains is what the typed Protein
                                       has. No flag, no marker.
```

The canonicalisation log (per `naming-canonical-vocabulary-2026-04-30.md` §7 schema) is the **methodology trail** — JSON record of what was read and what was applied, citable to upstream-authority files for reviewer spot-checking. The typed substrate is the **load-bearing state** that calculators read.

This is "nearly the only translating we do" because everything else is either (a) C-bucket label canonicalisation handled by the walk, or (b) data already in canonical form on disk that flows through unmodified.

**Morning task on this:** confirm this matches the user's earlier decision; refine if any specific calculator needs more than the existing fields offer; then NC1-NC7 implementation populates these existing fields from the upstream authorities.

---

## Water and ion access — walk-back of yesterday's over-abstraction

User flagged end-of-day:

> I made a bad choice yesterday and we talked about abstracting out waters and frame information. The truth is, this little tool should probably not be in the business of pulling out a bunch of frame data, waiting for the calculator running milliseconds later, and providing it back. The calculator can do ugly string stuff with waters and ions, they are only being looked at to get out electrostatics and exposure. But that does mean frame extraction has "how to get back to that frame" or some kind of cache that lasts for the whole frame pass.

**Decision (recovered / walked back to):**

```
SCOPE                           Waters and ions are accessed only for
                                electrostatics (per-water E-field / EFG
                                from charges) and exposure (hydration
                                shell geometry, half-shell asymmetry).
                                Narrow special-purpose use.

NO TYPED SUBSTRATE               Waters and ions are NOT in the typed
                                Protein / Residue / ConformationAtom
                                substrate. The protein-chemistry
                                canonicalisation discipline does NOT
                                extend to waters and ions.

CALCULATORS ACCESS RAW          "Ugly string stuff" allowed: water-
                                and ion-touching calculators
                                (WaterFieldResult, HydrationShellResult,
                                HydrationGeometryResult, indirectly
                                ApbsFieldResult via grid) read OW,
                                HW1, HW2, NA, CL etc. directly. No
                                in-tool translation to canonical
                                O/H1/H2 names.

ACCESS PATTERN                   Per-frame cache OR "current frame"
                                handle. Lives for the duration of the
                                calculator pass on a single frame.
                                When the trajectory rolls to the next
                                frame, the cache invalidates.

                                nmr-extract is NOT in the business of
                                pre-extracting all water positions to
                                a typed buffer for the model to consume
                                later. Extraction's job is per-frame
                                typed protein features; water access is
                                by-side-channel for the few calculators
                                that need it.

WHAT THIS DEPRECATES             naming-canonical-vocabulary-2026-04-30.md
                                §6 "Water and ion canonicalisation" — the
                                aggregate canonicalisation log entries
                                for waters/ions, the SOL/OW/HW1/HW2 →
                                HOH/O/H1/H2 translation tables, the
                                virtual-sites-for-future framing — all
                                overkill given the narrow use case.
                                Simplify §6 to a short "raw access via
                                per-frame cache" note.

                                NC7 in the implementation order shrinks
                                substantially (or merges into a separate
                                "frame access pattern" doc that names
                                where the cache lives and how calculators
                                use it).
```

**Morning tasks on this:**

1. Confirm the walk-back matches the user's intent (the recovered shape above).
2. Decide where the per-frame cache or current-frame handle lives architecturally — on Trajectory? On ProteinConformation? On a new FrameAccessor?
3. Review the existing water-touching calculators (WaterFieldResult, HydrationShellResult, HydrationGeometryResult, ApbsFieldResult) to confirm they already read raw or need a small access-pattern adjustment.
4. Update naming-canonical-vocabulary-2026-04-30.md §6 to the simplified shape; shrink NC7 in §"Implementation order"; add or extract a "frame access pattern" subsection naming the cache mechanism.

This walk-back simplifies a lot of work that was queued under the typed-substrate canonicalisation umbrella.

## Topology next steps — complete inventory

Consolidated for clarity. Every topology-related work item in the queue, in dependency order. This is "real physics topology knowledge" being landed step by step into the typed substrate.

```
1. AMBER charge slice (six steps GREEN, uncommitted)
   The charge-source plumbing the rest depends on. Land per-step
   commits per the central plan's discipline.
   Status: 62/62 tests passing; awaiting commit.

2. PHASE 0 / NC1-NC7 — trajectory loader naming canonicalisation
   Walk topol.top rtp + sources/structure_extracted PDB; populate
   existing typed fields; emit canonicalisation log.
   NC1   Define typed CanonicalizationRecord on ProteinBuildContext
   NC2   topol.top rtp comment parser (pure C++)
   NC3   Structural-match against canonical reference PDB
   NC4   Wire walk into FullSystemReader::BuildProtein
         (clears the "unknown residue HISH" error)
   NC5   Tests: 1P9J + 1Z9B walks produce correct logs
   NC6   Extend log to other load paths (--orca / --mutant / --pdb /
         --protonated-pdb / --amber-prepared)
   NC7   Water and ion canonicalisation — SHRUNK per today's walk-back;
         see "Water and ion access" section above. Becomes "frame
         access pattern" doc rather than canonicalisation discipline.

3. Real-chemistry capture (decided; see section above)
   Use existing typed fields. NO new fields.
     Residue.protonation_variant_index    (HID/HIE/HIP/CYX/ASH/GLH/LYN/...)
     Residue.terminal_state                (N-/C-term variant)
     CovalentTopology bond list filtered
       by BondCategory::Disulfide          (SG-SG inter-residue bonds)
     ConformationAtom.partial_charge       (varies per variant via lookup)
     atom-set                              (implicit; only present atoms exist)

4. IUPAC topology debts (3 items, from 2026-04-26 IUPAC topology landing)
   a. pro-R / pro-S verification           Geometric assignment exists but
                                            unverified; need to verify
                                            against IUPAC 1979/1998
                                            reference geometry.
   b. variant overrides                    Non-standard residues need the
                                            override pathway (currently
                                            standard residues only).
   c. DFT-compare completeness             MutationDeltaResult currently
                                            stores only delta_shielding total;
                                            needs WT original + mutant
                                            compared + diamagnetic +
                                            paramagnetic, all six tensors
                                            plus three deltas. Can adopt
                                            AtomReference for matching.

5. Gemmi adoption (agreed in principle 2026-04-26, NOT YET implemented)
   Phase 1   Read-only CCD validator. NPY output. No internal data flow.
             Used for fleet-curation spot-check; not in the trajectory
             load path. Validates per-protein topology against CCD;
             flags discrepancies.
   Phase 2   Non-standard residue scaffold. Extends typed substrate to
             handle modified residues, glycosylation, ligands when first
             encountered. Currently rejected at the AMBER charge slice
             resolver boundary (fail-loud); Phase 2 lifts that.

6. PHASE 1 substrate buildout (N1.A through N1.G)
   ONLY after PHASE 0 closes. The COMPLETE chemistry+geometry substrate
   built before any existing calculator is migrated.
   N1.A   Substrate parser pass (invariant; no library; typed atom-name
          enums; LegacyAmberTopology)
   N1.B   CCD ingest (invariant; gemmi; stereo_config, R/S/pro-R/pro-S,
          formal charges; LegacyAmberTopology)
   N1.C   pH / ionization layer (invariant for our MD; protonation_state_
          at_pH, net_charge_at_pH, isoelectric_point; LegacyAmberTopology)
   N1.D   Per-conformation geometric topology (calculator output; ring
          geometry, peptide ω, disulfide χss, rotamer classification;
          per-frame on ProteinConformation)
   N1.E   RDKit perception (invariant chemistry; SMARTS, functional groups;
          cross-validate N1.B; LegacyAmberTopology)
   N1.F   Projection surface (functions on LegacyAmberTopology — the
          crystal projection rule made concrete)
   N1.G   The substrate as MD-condition specification (the typed
          substrate IS what an MD run consumes; per-residue
          protonation_state_at_pH + net_charge_at_pH + isoelectric_point
          are the loadable state)

7. PHASE 2 calculator migration (after PHASE 1 closes)
   Calculator-by-calculator off the legacy string surface onto the
   typed substrate. Drift OK, surprise not. Per-calculator review.
   Calculators that don't need the richer substrate stay on legacy
   indefinitely; no forced migration.

8. OpenBabel exit (consequence of PHASE 2)
   Not a separate phase. Current OpenBabel-based CovalentTopology is
   replaced by typed-substrate-derived bond information as PHASE 2
   migrates calculators. The exit happens as the dependency graph
   clears.

9. N4 — Post-migration cleanup
   The legacy string surface fully retired. Anything depending on
   it is either migrated or explicitly archived.

DEPENDENCIES

  Slice (1) → PHASE 0 / NC1-NC7 (2) → PHASE 1 (6) → PHASE 2 (7) →
    OpenBabel exit (8) → N4 (9)

  Real-chemistry capture (3) sits inside PHASE 0 / NC1-NC7 — it's
    the typed-field-population shape the walk uses.

  Water/ion access (NC7) is now light per the walk-back; not a
    blocker for anything else.

  IUPAC debts (4) are independent; can land alongside or after
    PHASE 0. The pro-R/S verification + variant overrides are
    upstream of N1.B (gemmi CCD ingest); the DFT-compare
    completeness is independent.

  Gemmi Phase 1 + 2 (5) sit before / alongside N1.B and N1.E.
    Phase 1 is read-only spot-check; Phase 2 unlocks non-standard
    residue handling.
```

## Sample MD runs in flight overnight

Per `spec/plan/md-session-handoff-2026-04-30.md`:

```
Two sample proteins              Same two as the first batch (BMRB-PDB-extracted,
                                   with published S² / T1 / T2 NMR data)
Replica shape                    1 × 15 ns single replica per protein
Cadence                          40 ps in extraction (nstxout-compressed=10000
                                   + stride 2)
Force field                      ff14SB / TIP3P
Per-protein wallclock target     ~25 minutes for ~25 K-atom proteins on a 5090
                                   (per the prior fleet's md.log baseline)
Expected morning state           At least one run completed
```

These runs prove Gate 1 of two for the 685-fleet redispatch: MDP setup + extractor pipeline land cleanly. Gate 2 (OF3-direct structures clean) is separate.

Whether the runs landed cleanly tells the morning session whether the prep + MDP package is right. The actual H5 schema validation (compression, drop aim, float32 positions) is downstream of the runs and a separate (not-prep) session's work.

---

## OpenAI critiques deferred for morning review

User said: "there were a couple of openai critiques we put off for good reasons, and then hit it with another review."

Likely the review files in `spec/plan/`:

```
spec/plan/openai-5.5-strong-architecture-layout.md
spec/plan/openai-round-4-response-2026-04-28.md
spec/plan/opus-round-2-review-2026-04-28.md
spec/plan/opus-round-3-review-2026-04-28.md
```

User remembers what was put off and why. Recover the specific critique items by reading these files; don't try to reconstruct from session memory.

---

## Suggested order for the morning

1. **Read this file first.**
2. **Check whether sample MD runs landed cleanly overnight.** If yes, Gate 1 passed; the prep + MDP package is validated.
3. **Review the deferred OpenAI critiques** in `spec/plan/openai-*.md` and `spec/plan/opus-*.md`. Identify which items the user wants to revisit now.
4. **Refine the topology-capture question** above. Land typed substrate fields on Residue / ConformationAtom for protonation variant + disulfide partner + atom-added / atom-stripped flags.
5. **Bridge the docs.** Update `amber-implementation-plan-2026-04-29.md` working-tree-status block to reflect today's PHASE 0 partial + naming spec + open topology question. Update `naming-canonical-vocabulary-2026-04-30.md` to reference the central plan's commit-ready state. Update `spec/INDEX.md`, `KNOWN_BUGS.md`, `CLAUDE.md` reading-order as needed.
6. **Commit.** Per-step commit discipline for the slice (six commits per the central plan); separate commits for the naming spec, the PHASE 0 partial, today's planning docs, and the spec-doc revisions.
7. **Start NC1.** Define the typed `CanonicalizationRecord` on the substrate; this is the first actual code step of NC1-NC7. NC4 is the step that clears the "unknown residue HISH" error.

---

## Key docs by purpose (the morning's reading)

```
ABOUT WHERE WE ARE
  spec/plan/session-handoff-20260430-evening.md      this file (start here)
  spec/plan/running_plan_notes.md                    today's running record

THE SLICE + PHASE 0 + NAMING WANDER
  spec/plan/amber-implementation-plan-2026-04-29.md  central plan
                                                       (working-tree-status
                                                       block stale; will be
                                                       updated in step 5)
  spec/plan/naming-canonical-vocabulary-2026-04-30.md naming spec (NC1-NC7
                                                       implementation order)
  spec/plan/current-topology-anchor-2026-04-29.md    topology anchor
  spec/plan/pre-iupac-cruft-map-2026-04-29.md        pre-iupac cleanup map

THE MSC FRAMING
  spec/plan/thesis-discussion-part1-2026-04-30.md     MSc framing
  spec/plan/comprehensive-calculator-inventory-...md  calculator inventory
                                                       (vetted, advisor-shareable)

MD WORK
  spec/plan/md-rerun-685-discussion-priors-...md      MD priors
  spec/plan/md-session-handoff-2026-04-30.md          prep handoff
                                                       (the brief that's
                                                       producing the runs
                                                       overnight)

DEFERRED REVIEWS
  spec/plan/openai-5.5-strong-architecture-layout.md
  spec/plan/openai-round-4-response-2026-04-28.md
  spec/plan/opus-round-2-review-2026-04-28.md
  spec/plan/opus-round-3-review-2026-04-28.md
```

---

## What this file is NOT

Not a replacement for any of the above docs. It is a pickup brief — read it first, then go to the docs it points at.

Not a plan-of-record for any technical decision. The naming-canonical-vocabulary spec is the plan of record for the trajectory loader naming work; this file just notes that it exists and where it sits in the sequencing.
