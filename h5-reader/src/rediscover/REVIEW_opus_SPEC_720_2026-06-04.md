# Opus adversarial review — SPEC_720_STATIC_INGEST_2026-06-04.md

> **Historical review record — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Reviewer: opus, adversarial-but-fair. Scope: logic + consistency of the spec
against the repo + the source brief. Docs-only; no code/git/runs. Read-only
everywhere else. Branch `h5-reader-pysr-spike` — lead owns all git.

Verdict up top: **the spec is strong and faithful to the brief.** Grounding
citations are real and honest (I spot-checked the load-bearing ones below).
The foreclosure is essentially airtight, with two small soft spots to harden.
The memory-strategy fork is kept honestly open and correctly tied to the
unsettled between-axis maths. Findings follow in brief order.

---

## 0. Citation audit (did the spec tell the truth about the code?)

Verified against source, all confirmed:

- `SingleConformation` is real, reads `Pos` from the snapshot, `frameCount()==1`,
  `timePicoseconds()==0.0`, explicitly "NOT a faked one-frame trajectory"
  (`../model/SingleConformation.h:1-43`, `.cpp:18-30`). Spec §2 accurate.
- `QtProteinLoader::LoadPose` exists for `SinglePose`, loads sidecar +
  flat per-atom NPYs via `FrameNpyLoader`, builds `SingleConformation`
  (`../io/QtProteinLoader.h:64-86`). Spec §2 accurate.
- `FrameNpyLoader` is lenient exactly as described: enumerates `*.npy`,
  skips unrecognized/malformed, returns partial snapshot, null only on a
  dir-level failure (`../io/FrameNpyLoader.h:17-21`, `.cpp:46-122`). Spec §2/§13
  accurate — **and see §1 finding below, this is the most important seam.**
- `CalcsetManifest::ResolveLgsPath` resolves a dir by listing a single `*.LGS`,
  hard-errors on 0 or >1 (`../io/CalcsetManifest.cpp:147-178`). Spec §4's
  recommendation to pass exact `.LGS` paths (avoiding even this single-match
  enumeration) is sound.
- DFT loader 0.1 ppm identity tolerance is real (`../io/DftShieldingLoader.cpp:39`,
  `kIdentityTolPpm = 0.1`, checked on T0 at `:104`). Spec §11 gate cite accurate.
  NOTE: the identity is validated on **T0 only** ("decomposition is linear, T0
  stands in") — this is exactly the maths-audit issue #1 in `NOW.md:49`. Fine for
  the gate, but see §4 finding.
- `OrcaShieldingParser` stores `total_raw/dia_raw/para_raw` 3×3 + `orca_coord`
  (`../io/OrcaShieldingParser.cpp:99-122`). Spec §5 target table accurate.
- `PerAtomSubstrate` size-gate is real and throws (`PerAtomSubstrate.cpp:3786-3790`).
  Spec §10 reuse claim accurate — **and the spec sharpens it** (gates on TOTAL
  output, not just the append slab; see §5).
- `NODE_STORE_CONTRACT_2026-06-02.md:35` states "THE LAW: no maths model 2 in
  Python." `PARTITION_FILTER_DESIGN.md:16-23` states the C++/Python boundary
  exactly as the spec quotes. Materialization-as-named-transient is real
  (`NODE_STORE_CONTRACT:96-101`). Spec §1 grounding accurate.

**Line-ref nits (cosmetic, fix if cheap):**
- §3 / §13 cite the trajectory gate at `RunData.cpp:84`. The actual reject is
  `RunData.cpp:86` (`manifest.kind != Trajectory`). Off by 2; the `:76`/`:87`
  brackets elsewhere are fine.
- Several `../io/QtFieldCatalog.gen.h:NNN` line cites in the §5 table were not
  re-verified line-exact (the field groups exist; only the exact lines are
  unchecked). Low-stakes.

---

## 1. THE FORECLOSURE — is it airtight? (the #1 check)

**Largely YES.** §1 names both failure modes as LAW 1 / LAW 2, not advice; the
"Rationalization guard" pre-empts the "I had to do it this way" move and routes
it to a C++ extension on **both** memory branches (§1, §7). §7's recommended
structure keeps proteins/indexes in C++ for branch (a) and only C++ reduced
buffers for branch (b) — neither hands Python a model. §12's Python gate is
explicit: fitter consumes only `static_rows.*` / `static_*.npy` / specs /
manifest, may not import reader model code, open `.LGS`/`trajectory.h5`/ORCA,
or scan producer dirs; frozen `get_C` honored. This is the right shape and it
structurally forecloses both the terabyte and the second model.

**Two soft spots to harden (neither fatal, both worth a sentence of LAW):**

**(1a) The glob seam is the one real leak, and it's under-nailed.** The spec's
own reuse target, `FrameNpyLoader::LoadSnapshotDir`, does
`d.entryList("*.npy")` at `../io/FrameNpyLoader.cpp:46` — that **is** file
discovery (enumerate-and-match), the exact thing `feedback_no_file_discovery`
forbids, and it silently tolerates a missing required field (column just stays
absent, partial snapshot is "success"). §2 and §13 correctly say the static
path must NOT reuse this leniency and must use a "strict expected-field loader,"
and §4 lists the required per-protein checks. **But the spec never states as LAW
that the static loader resolves each required NPY by its documented exact path
and never enumerates the directory.** As written, a future builder could
"reuse `FrameNpyLoader` in strict mode" and keep the glob. Sharpen: add to §4 a
one-line law — *"the static loader constructs each required NPY path from the
documented field→filename convention and opens it directly; it never calls
`entryList`/glob over the pose dir. A missing required path is log-and-stop."*
This closes the only concrete no-file-discovery hole in the design.

**(1b) "Optional source-level / pair diagnostics" is a small crack the brief
would want explicitly sized.** §8 ("Drop As Trajectory-only") and §10
("FORECLOSED") both allow a *named, capped, C++-produced* pair/source diagnostic
"for an explicit gate," never the fitter substrate — consistent with
`NODE_STORE_CONTRACT:96-101`. That is the right carve-out and I would NOT remove
it. But §10's "even a temporary Python-readable pairwise dump for fitting" is
forbidden, while §8 allows a "named-query diagnostic." The boundary between
"capped C++ diagnostic" and "the 57 GB anti-pattern" is the cap. Sharpen: state
the cap as a number or a rule (e.g. *bounded by top-k sources × a small atom
subset, never all-pairs × all-atoms, and never an input to the fitter*), so
"named and capped" can't quietly grow. The existing `PerAtomSubstrateConfig`
already models this discipline (`reader_pair_atom`, `reader_pair_frame_slots=3`,
`top_k=3` — `PerAtomSubstrate.h:33-39`); cite it as the template.

Net: **foreclosure is airtight on the two headline modes; the residual risk is
the glob seam (1a), which the spec discourages but does not outlaw.**

---

## 2. MEMORY-STRATEGY FORK — honestly open + maths-tied?

**YES, this is done well.** §7 lays out (a) all-resident and (b) sequential
buffers, refuses to pick, and ties the choice explicitly to the unsettled
between-axis statistics method (§7 intro, §13 "OPEN - memory strategy",
"OPEN - between-axis statistics method"). The recommended non-foreclosing
structure — one `StaticCohortAccumulator` API with two backends, identical
`StaticRowDigest` in and identical `StaticEmitWriter` out, default unset until
the lead chooses — is exactly the right move: it forecloses neither branch and
keeps spatial+model in C++ on both. This is consistent with `NOW.md:88` /
`STATE.md:65-72` (between defers entirely to 720-WT; 1P9J has no clean between
axis), so the fork is correctly motivated, not decorative.

**One honest-openness sharpening:** §7 branch (b) says "If the maths needs all
row vectors co-resident, branch (b) may not be valid." Good. But the spec leans
faintly toward (a) by listing (a) first and richer ("maximal flexibility"). Keep
it symmetric: note that branch (a)'s cost isn't only RAM — an **all-resident
cohort statistic (global centering, cohort quantiles) is itself a between-axis
maths commitment** the lead has NOT made, and choosing (a) "to be safe" would
silently bias toward population-level reductions the between-axis discussion may
reject. i.e. (a) is not the free/safe default; it presumes a class of estimator.
That keeps the fork genuinely open rather than (a)-leaning.

---

## 3. CONSUME-THE-GIVEN path — concrete? no re-implementation creep?

**YES, concrete and disciplined.** §5's source table maps every payload to a
real producer surface + a real reader loader: identity←`.LGS`, topology←
`QtTopologySidecar` (5 NPYs), positions←`FieldKind::Pos`, ring/bond/McConnell/
electrostatic shadows←the producer NPY field groups, DFT target←raw ORCA
dia/para/total via `DftShieldingLoader`/`OrcaShieldingParser`. §6's per-protein
recipe (load `.LGS` → topology → strict snapshot → one-frame conformation →
target → catalog → `ResidentIndexes` → traverse) is the trajectory load shape
minus the time axis, which is precisely the brief's ask. The "no re-implement
nmr_extract, no kernel recompute reader-side" rule is stated as RESOLVED (§5,
§13) and the only C++-side computation is neighbourhood **geometry/reduction**
over producer positions — which is the model-is-spine contract, not a recompute
of physics. No creep.

**Two findings:**

**(3a) The producer-rerun GIVEN is doing heavy lifting and disagrees with an
older in-tree pointer — the spec is faithful to its brief, flag for the lead.**
The brief (`CODEX_BRIEF...:30-42`) makes "extractor re-run over all 720 curated,
MOPAC ON, producer output exists per protein" a GIVEN, and declares the 725 in
`/shared/2026Thesis/consolidated/` out-of-date (duplicates). The spec follows
this correctly. BUT `NEXT_SESSION_PROMPT.md:45` still says the 720 statics use
"the EXISTING 720 DFTs ... `/shared/2026Thesis/consolidated/`, NO new ORCA." So
there is a live contradiction between the brief (re-run, MOPAC ON, consolidated
is stale) and the standing continuation prompt (reuse existing, consolidated is
the source). The spec is right to obey the newer brief and stay path-agnostic
(§0, §4), but **the lead should resolve which is canonical before the build
loop** — the whole MOPAC-shadow column set (§5, §8) is contingent on the re-run
actually happening with MOPAC ON. The spec flags MOPAC presence (§4, §8 "mopac_on
flag"), which is the right hedge; just surface the brief-vs-prompt conflict
explicitly so it isn't silently inherited.

**(3b) Static target source is correctly left as a flagged fork, but note the
ORCA-coord parity dependency.** §5/§13 recommend raw ORCA via `DftShieldingLoader`
as canonical target, NPY arrays as cross-check. Good. The raw-tensor T2 is only
frame-comparable if the ORCA-input coords match the producer positions — the
existing tree already worries this (`ExtractionSupport.cpp` `CheckDftFrameAlignment`,
the Kabsch diagnostic; `DESIGN.md:80-83` "cross-frame T2 comparison unverified").
For statics this is even more pointed: there is no trajectory of frames to
average the rotation check over. T0 is safe (rotation-invariant). The spec should
note that the static **T2 target inherits the same frame-alignment caveat**, and
the per-protein manifest (§9) should carry the single-pose Kabsch angle so a
cross-frame T2 isn't trusted blind. (T0/|T2| safe regardless.)

---

## 4. ORACLE-PARITY GATE — sound + testable?

**YES, the strongest section.** §11 correctly separates Must-Match (exact
integer equality on identity/topology/source IDs/counts; tiered float
tolerances by origin — 1e-9 Å positions, 1e-12/bit-identical for copied f64
kernels, 5e-7 for f32-origin embeddings, 1e-9 ppm for same-text DFT vs 1e-4 for
NPY-rounded) from May-Differ (`row_id`, `h5_row`/`frame_slot`/`time_ps`,
trajectory-only Welford columns, run-kind/provenance). The join on
`(protein_id, atom_index)` ignoring row order is correct since static uses
cumulative offsets. The comparator-is-Python-only-over-emitted-rows constraint
(no opening `trajectory.h5`/ORCA/sidecars/producer NPYs) is the right teeth and
is consistent with `PARTITION_FILTER_DESIGN.md:16-23`. The tiered tolerances are
physically motivated (f32 producer arrays → 5e-7), not magic numbers — meets
`feedback_transparent_cutoffs`.

**Two sharpenings:**

**(4a) Define what "the same coordinates" means, precisely.** §11 Gate Input
says "a static-pose producer output generated from the same coordinates/topology
and same raw ORCA target" as one 1P9J trajectory frame. The gate's whole
validity rests on the static pose being **bit-identical coordinates** to the
chosen trajectory frame, not merely "close." State it: the oracle static pose is
produced from the exact extracted positions of a named trajectory frame, so
positions must match to f64 round-trip (`abs ≤ 1e-12` / bit-identical), and the
1e-9 Å budget is reserved for derived geometry, not the input positions. As
written, "abs ≤ 1e-9 Å" on positions (§11) is looser than the gate can afford if
the point is to prove the static path reproduces the trajectory path on
identical input.

**(4b) Scope the gate to what `--case all` actually covers, or say it doesn't.**
`NODE_STORE_CONTRACT:145-147` records that the existing oracle-parity command
covers **ring + mc parity only** — NOT broad/all-atom/efg/buckingham/aimnet2.
§11's feature-family parity list (§11 "Feature-family parity") is much broader
(ring, charge, McConnell, MOPAC, APBS, AIMNet2, hbond/pq/disp/HM/ringchi/water/
SASA/EEQ). That's the right *aspiration*, but the spec should note the existing
comparator does not yet cover those families, so the static oracle gate either
(i) extends the comparator C++-side to emit per-family max-diff for the full
menu, or (ii) explicitly states which families are gated v1 vs deferred. Right
now §11 implies full-menu parity that the current tooling can't deliver — an
overclaim by omission. Cheap fix: one caveat sentence + a pointer to extend the
emitted parity report, not the Python comparator.

---

## 5. CONSTRAINTS — model-is-spine / <15 G / extractor SACRED / no file
discovery / frozen `get_C`?

- **model-is-spine:** honored. §1 LAW 2, §6 (all neighbourhoods are C++ verb/
  index queries), §7 (both branches C++-side), §12 Python gate. Consistent with
  `GUIDANCE.md:22-34`, `feedback_model_is_spine`. ✅
- **lean disk < 15 G:** honored and well-argued. §10's arithmetic is explicit
  (R × cols × bytes), spans A_avg 1k/2k/4k (1.8 / 3.6 / 7.1 GB without embedding;
  worst-case 10.1 GB with), and the preflight size-gate refusing >15 GB reuses a
  real fail-loud pattern (`PerAtomSubstrate.cpp:3786`). **Sharpening:** the
  existing gate fires at 10 GiB on a column *subset*; the spec proposes 15 GB on
  *total*. That's a genuine improvement, but flag the mismatch so the builder
  doesn't copy the 10 GiB-subset gate and think it's done — the static gate must
  sum ALL enabled sidecars + identity + manifests. Also: §10's C64 ≤ 220 default
  is a lot of columns for a one-frame substrate; the column-count constants in
  `PerAtomSubstrate.h:131-141` (classical 89 + conditioning 32 + dominance 10 +
  ring-path 226 + method-path 111 + hbond 73 …) sum well past 220, so "carry
  these families" (§8) and "C64 ≤ 220" (§10) are in mild tension. Not wrong —
  §8 is explicitly lead-vettable and smaller-than-trajectory — but name the
  tension: the static menu must be a genuine *subset*, and 220 is a budget
  ceiling, not a carry-everything license. ✅ with the caveat surfaced.
- **extractor SACRED:** honored. §0/§5/§13 — consumes producer output, does not
  re-run or modify `nmr_extract`. Consistent with `GUIDANCE.md:17-20`,
  `h5-reader/CLAUDE.md` ("never link the library, never trigger extraction").
  One note: the static ingest lives in the **reader/rediscover** tree and reads
  the producer's NPY surface — it does not touch the library extractor. The spec
  is clean here, but the wording "re-run the extractor over all 720" (the GIVEN)
  is the *producer's* action, upstream and external; the rediscover ingest never
  invokes it. The spec already keeps this straight; worth a half-sentence so no
  reader thinks rediscover triggers extraction. ✅
- **no file discovery:** honored in *intent* (§4, §13) but see finding **1a** —
  the glob in the reuse target is the one place the design discourages-but-
  doesn't-outlaw enumeration. Harden §4 to LAW. ⚠️→fixable.
- **frozen `get_C`:** honored. §12 Python gate: "Frozen `get_C` remains the basis
  bridge on the Python/e3nn side; no Python recomputation of DFT tensor
  decomposition or source geometry." Consistent with `NODE_STORE_CONTRACT:35-38`.
  ✅

---

## 6. Are the flagged forks the RIGHT ones, and honestly open?

**Yes.** §13's Open/Caveated list captures the real decisions: memory strategy
(tied to maths), final static schema (embedding/T1/solvation inclusion), canonical
DFT target source, exact cohort manifest shape, between-axis statistics method,
output compression, and the two CAVEATS (exact disk/RAM counts ungrounded by
design since path is TBD; `FrameNpyLoader` leniency must not be reused). These
are the genuine unknowns, and none is prematurely resolved.

**Missing-from-the-list forks the lead will want surfaced (add to §13):**
1. **Brief-vs-prompt conflict on the input source** (finding 3a) — re-run-with-
   MOPAC vs reuse-existing-720-DFT. This is upstream of "exact cohort manifest
   shape" and bigger than a path fill: it decides whether the MOPAC columns exist.
2. **Static T2 frame-alignment caveat** (finding 3b) — does the single-pose target
   carry a Kabsch angle, and is cross-frame T2 trusted? T0 is safe; T2 is not free.
3. **Oracle-gate family coverage** (finding 4b) — full-menu parity vs v1 ring+mc,
   given the current comparator only covers ring+mc.

---

## Bottom line

The spec is honest (BUILT/DESIGNED/OPEN/CAVEATED is used with discipline and the
code cites are real), faithful to the brief, and gets the two hard things right:
the foreclosure is structural (LAW 1/LAW 2 + the rationalization guard + the
Python gate), and the memory fork is genuinely open and maths-tied. The residual
work is hardening, not redesign:

1. **Outlaw the glob (1a)** — make "resolve each required NPY by exact documented
   path, never enumerate the pose dir" a LAW in §4. This is the single most
   valuable change; it's the only concrete no-file-discovery hole.
2. **Surface the brief-vs-prompt input conflict (3a)** and the static-T2 frame
   caveat (3b) — both are silently-inheritable and both gate real columns.
3. **Right-size the oracle gate (4b) and the column budget (§5)** — say the
   comparator covers ring+mc today and the full menu is a C++ parity-report
   extension; name that C64≤220 is a ceiling, not a carry-everything license.

None of these blocks the lead's plan-vet; they're the sentences that keep a
future build loop from re-deriving the boundary cold.
