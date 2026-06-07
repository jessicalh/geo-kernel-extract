# Continuity — kernel redesign & the run from data to model (2026-06-06)

*Draft for Jessica to shape. This is the session handoff across contexts: where the work
stands, the principle that governs it, what is still open, and the failure modes the next
instance must be warned about. Jessica takes the held-open issues next; they are flagged
here, not resolved.*

---

## The arc (the framing)

Our initial geometric-kernels conception was an interesting **experiment** — worth doing for
its simplicity and its back-of-the-napkin estimation power, the chance that a pile of
hand-computable forms might carry real signal. It did its job as an experiment. But now we
have the **framework, the data, and the toolchain** to run the whole way from data → statistics
→ model — and so we are **replacing those first kernels with better work**: properly grounded,
cited, equivariant feature designs (the `kernel_design/*.md` set), fed to e3nn, swept and
ablated under an honest stats program. The napkin estimate was the start of the story, not the
claim.

## The principle that governs inclusion (the underlying thing to remember)

Every feature may have a **scalar (`0e`) form** and a **tensor (`1o`/`2e`) form**. The rule:

> **A feature's scalar form, if it has one, can drop into the final "fido & biscuits" shift
> predictor if it makes sense there. A feature should go away only if it makes sense in
> NEITHER context.**

The two contexts are (1) the **shielding-tensor work** — the explanation/attribution side, where
the rank-2 content and the partition between kernels matter; and (2) the **final shift predictor**
— the big wild-west model against experimental shift, where overlap does not matter and a useful
scalar earns its seat regardless of attribution. So a feature that is confounded or redundant for
*attribution* can still ride in as a *scalar* in the shift predictor if it helps there. It is cut
only when it earns nothing in either place. This is the inclusion test, not "does it correlate
with DFT" alone.

---

## State of the work (pointers, not re-derivation)

- **Kernel-design program — complete for the core set.** Seven design docs in this directory.
  Honest tally: a **solid core of three** — ring, McConnell (bond anisotropy folded in, same
  Δχ⊗propagator→`2e` machine; extend `McConnellResult`), charge/EFG — plus three **carried-along
  maybes** the partition work decides: H-bond/Larsen (drop-candidate — the Larsen term is circular
  as a feature and its electrostatic overlaps charge/EFG; *carried, not cut*), π-quadrupole (may be
  charge/EFG's next multipole order of the same charges), dispersion (may not earn a slot). The set
  converges toward ~5; the arc is robust to the count; **cut by the two-context inclusion test, not
  a round number.** House pattern: clean equivariant spherical-tensor features (relevant
  `0e`/`1o`/`2e`), raw in the molecular frame, geometry unscaled, scales fitted or as parallel
  channels, e3nn forms the invariants, ceilings disclosed.
  Docs: `ring.md`, `mcconnell.md`, `charge_efg.md`, `bond_anisotropy.md`, `dispersion.md`,
  `hbond_larsen.md`, `pi_quadrupole.md`.
- **MOPAC — reclassified to first-class, run on everything (per-frame).** No longer the
  retired/static FullFat probe. Because the OpenFold structures are poor enough that MD must be run
  on everything anyway just to earn ~1 ns of good MolProbity, MOPAC rides along at marginal cost and
  **will be run per-frame on everything**; its output was *strong* in testing. It is a first-class
  **local-electronic** source beside the through-space geometric kernels: conformation-responsive
  Mulliken charges (a strong charge channel), Wiberg bond orders (the McConnell weighting), and s/p
  orbital populations + dipole (a local-electronic descriptor class the geometric kernels structurally
  cannot reach). **Full-capture shape now spec'd** (`mopac_extension.md`): we'd been
  emitting a *sliver* — four arrays from a fraction of the `.out` — of a run that prints vastly more
  (full MULLIK incl. AO/shell populations + overlap populations, the full directed bond-order object
  with valency diagonals, the full dipole table, the energy decomposition, MOZYME localized orbitals,
  plus the raw `.out`/`.aux` as lossless backstop). **Implementation in flight** (extend `MopacResult`
  through the SDK; full capture held in memory — 128 GB, generous RAM, lean disk; deprecate-and-add
  over the four legacy arrays; respect ConformationAtom + Result conventions; NOT the reader, its own
  effort). **Sidecar format is discuss-first** (codex proposes, human reviews). Partition caveats:
  MOPAC charges overlap AIMNet2's; MOPAC-bond-order McConnell sits beside literature-Δχ McConnell —
  channels the fit weighs, overlapping for attribution.
- **Pipeline adaptation** (`pipeline_adaptation.md`): the substrate already carries irrep fields
  (`ArrayRank` Vec3/T2_5/Tensor9; `value`/`valueVec3`/`valueT2`/`valueTensor`), so the change is
  mostly **registration** plus the attacher/schema in `ComposedRelationships.cpp`. The geometry
  layer is good and untouched. **Naming is the load-bearing part and is decided with Jessica, not
  autogenerated** — name-is-the-thing for both the stats pass and the functional interface.
- **Stats program & guardrail** (`doc/emerging/STATS_PROGRAM.md`): say the sensible things at the right
  confidence and leave the uncertain ones alone; the paper is not the code; judge against
  *sensible*, not *proven*; do not load the dice (the analysis hopes nothing); per-element AND
  per-IUPAC charts every time, even at small n; the sweep = complete PySR → equivariance where it
  makes sense → always an R², run-all-then-select.
- **What the kernels are** (`doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md`): low-order (l≤2)
  multipole moments of the weighted neighbourhood; physical-field and geometric-descriptor are the
  same object; the **circularity** — feeding the gradient field in makes form-recovery tautological,
  only the DFT-target match is real.
- **The open description** (`doc/emerging/WHAT_ARE_WE_DOING.md`) and the **lineage / borrowed
  theorems** (`doc/emerging/GEOMETRIC_KERNEL_MATH_LINEAGE.md`, `doc/emerging/STAIRCASE_SOCIAL_HISTORY.md`).

## The rewrite architecture (the plan)

The integration is **one new encompassing calculator** — honestly named *MOPAC, extended*. It is
added **beside** the existing calculators, which **stay in place as the exploration** (the napkin
experiment is the record of how we got here, not deleted). Producer discipline: **vet that the new
calculator leaves the old ones untouched** — they keep emitting exactly as before.

Shape:
- **Foundation** = everything the geometry layer already knows (neighbourhoods, frames,
  displacements — the good layer, kept) + everything MOPAC knows (charges, bond orders, s/p
  populations, dipole). **MOPAC is treated as ground truth** — the authoritative local-electronic
  information. The calculator is, honestly, *an extension to MOPAC plus our geometry*, not "our
  shielding model."
- **Composition** = every feature (ring, McConnell, EFG, …) built on that one shared substrate as a
  named composition, so shared sources are visible and each partition (e.g. π-quadrupole vs
  charge/EFG over the same ring charges) is decided **once at the composition**, not chased across
  files.
- **Default-source rule.** The new calculators **default to MOPAC** — we are a MOPAC extension — and
  **override with ours only where ours is demonstrably better.** The obvious overrides aren't close
  calls: **geometry is ours** (MOPAC doesn't frame positions/neighbourhoods/frames) and **the
  through-space kernels are ours** (MOPAC doesn't compute a ring current); MOPAC owns the
  local-electronic (charges, bond orders, populations). So the split is mostly clean by construction —
  MOPAC for local-electronic, us for geometry-and-through-space. For genuinely-overlapping quantities
  (a charge both have), MOPAC is the default and the override must be **shown, not assumed.** This is
  also the credible stance: *extend a real QM method, override only where demonstrably better* beats
  *our home-rolled values are the truth.* The overlap inventory flags the choices; this rule resolves
  them.
- **God object inside, normal calculator outside.** Internally it holds the orbicular view of the
  whole protein and all desired outputs — the point, the thing the narrow old calculators cannot do.
  Externally it is just a Result + a set of NPYs like the others; the god-object-ness does not leak
  into the interface.

Spec items:
- **MOPAC runs first** — the new calculator depends on it, so the ordering is a real dependency in
  the spec, not an assumption.
- **Leave the old calculators alone** — vetted; the old keep their exploratory emit unchanged.
- **Propagate MOPAC to the reader model** — the useful per-bond and per-atom MOPAC information (bond
  labelling/orders, charges, orbital populations) propagates out to the reader's model
  (`QtAtom`/`QtBond`) with this change, so the UI can show **what MOPAC knows on the model when an
  atom or bond is selected.** The reader is the vetting surface; and because MOPAC is now per-frame,
  those values are also time-series the illustrator can carry over the trajectory. Part of the plan,
  not polish.

### Held-whole vs context-sized — the resolution

The tension: calculators were context-sized for a reason — a whole held in one context sees the
connections, and it shows in the code — but the MOPAC synergy must NOT be fragmented across pieces.
Resolution: **hold MOPAC whole in a buffer; fan context-sized agents around it.**
- A **MOPAC buffer** ("MOPAC Protein") snarfs everything MOPAC knows, squirrels it away, and
  **announces it to each agent coming in** — held whole, undivided, in one typed place.
  Buffers-from-the-ctor pattern: attached at protein/trajectory scope, sized at construction,
  per-frame facts (MOPAC runs on everything now). Honours the Protein wall (a separate buffer, not
  geometry on Protein); objects-answer-questions made literal. **Its interface names are
  load-bearing** — every per-metric agent reads MOPAC through them, so they're decided deliberately.
- **One agent per metric** (T0/T1/T2, NPYs out, as always) reads the whole buffer instead of
  re-deriving — context-sized, with the full MOPAC view available. **Worktree them** for parallel
  isolation.
- Coordination runs through the **spine (the buffer)**, not through documents/process — the kind
  that shows in the code, the kind humans don't fail at.

**Implementation decision (2026-06-06, confirm reading):** do NOT build a brand-new MOPAC-Protein
object + sidecar by EOD. **Extend the existing `MopacResult` into the held-whole buffer role and
document the heck out of it** — it already snarfs MOPAC; extend it to hold/serve everything, heavily
documented (the "comes with the nice docs" point).

**The hard part — retiring the old MOPAC path, and the one exception to "old stays as exploration."**
MOPAC runs **once** (expensive; no cheap re-run to fix a mistake). Everywhere else old and new coexist
because re-running is cheap, but the old scattered MOPAC consumers (`MopacCoulomb`, `MopacMcConnell`,
the reconciliation probe) cannot keep drawing from the single expensive run *and* sit beside the new
path without the wiring fighting. The run-once economics **force consolidation, not parallel-keep** —
so retiring/rewiring the old MOPAC consumers carefully (without re-running MOPAC) is the real work, and
the place the careful spec passes earn their keep.

## Held-open issues (Jessica takes these next — flagged, not resolved)

- **The partition web.** The kernels overlap; almost every one's honest uncertainty is now a
  partition decision against another (bond-anisotropy↔McConnell, aromatic-bonds↔ring,
  dispersion↔ring, H-bond↔charge/EFG, π-quadrupole↔charge/EFG and ↔ring, plus MOPAC-charge↔AIMNet2
  and MOPAC-bond-order-McConnell↔literature-Δχ). These partitions are load-bearing and mostly
  **unsettled in the code**. Overlap is fine for the predictor, a confound for attribution. Do not
  resolve by fiat.
  - **An inventory we need (queued follow-up, can wait), NOT the partition evidence:** a separate
    agent, run *after* the MOPAC-extension capture comes back and shares — *what does MOPAC tell us
    that our existing geometry and work tell us differently?* The full map of MOPAC↔geometry overlaps,
    once the complete MOPAC output is in hand. **Not to exclude any MOPAC output — to know the
    overlaps, so we reach for the right source when composing the new calculations.** Be precise: this
    is *not* the partition evidence. The partitions are properties of the **new** calculations, which
    change (they use MOPAC differently), so they cannot be read off the existing overlaps. The
    inventory is **upstream of the build** — which source feeds which calculation; the partition is
    **downstream**, a property of the built new calculations.
- **What the funky kernels encode** is genuinely open (clean for the l=2 ones, open for the
  richer/empirical ones). "Open" is the true word.
- **Dispersion may not earn a slot** — no separable dispersion shielding to anchor against; a
  marginal-fit `0e` surrogate at best. The sweep decides, not us.
- **Why MOPAC tested strong (held open — plausible account, not asserted):** the geometric kernels
  are through-space and l≤2, structurally blind where shielding is local-electronic-dominated (the
  heavy-atom / "H-was-the-only-fair-test" case); MOPAC's electronic descriptors are a cheap proxy
  for exactly that local structure. So MOPAC may be the *local complement* carrying what the
  geometric set leaves on the table. Test: whether MOPAC's lift concentrates on the heavy,
  local-dominated atoms. Do not claim it; verify it.
- **Retractions/recharacterizations to carry, not re-litigate:** 1P9J between-atom LOAO is
  ~null → 1P9J is a *within-instrument* result, transferability deferred to the **720-WT** corpus;
  the EFG→DFT-T2 echo collapsed under a lab→local frame correction (it was a rotation confound).
  Frame as scope clarified, not as loss.

## Process state

- **Not ready to integrate.** When we do, it is multiple codex agents, worktrees, and one solid
  all-at-once hit — *because it must be* (the emit, reader, and e3nn input change together). The
  honest accounting here is the prerequisite; the hit can't be taken twice.
- **This is becoming a rewrite, not incremental deprecate-and-add.** The kernel set is a new,
  coherent design (proper-lego + MOPAC-anchored), not the old set patched; provenance lives in git
  and the deprecated design docs, not in dead kernels kept live in the emit. **MOPAC, run on
  everything, is the reliable anchor the rest hangs off** — it moves from a conditional
  "absent-not-faked" source to always-present, so kernels and descriptors may depend on it.
  **Surgery (clean names in place) still governs the reader/functional interface** — name-is-the-thing,
  accumulation corrupts there.
- **Fixed-DFT constraint:** the DFT campaign is done and not re-runnable (six months; shielding
  only — no magnetizability/susceptibility). Designs use geometry + literature/tabulated constants
  (Δχ from tables, fitted) + what the DFT already wrote. Re-extraction is cheap (<1 day); the DFT
  is the precious fixed thing.
- **References:** foundational primaries (Haigh–Mallion 1972, McConnell 1957) are frayed out — not
  openly posted anywhere legitimate. **H&M is covered by held secondaries** (Case 1995, Moyna 1998,
  Gomes–Mallion 2001); **McConnell's strong secondary** is the Case review at PMC3877577 (browser
  pull — PMC has a proof-of-work wall). The kernel-design citations were grounded on search
  summaries (WebFetch was limited) — **verify each on the pull** before any is load-bearing in prose.

## Failure modes of the tool (co-named — carry forward and correct)

These are where the AI is least trustworthy on this work; the protection is structural, not trust:
the honest description lives **on the page** (survives context loss), claims are **verified from
code** not asserted from memory, and Jessica reviews hardest at these zones.

- **Under-rates the faint and the funky** — calls "sensible with a ceiling" a "null, drop it."
- **Collapses open ambiguity** into a clean label (descriptor-vs-mechanism, the round-trip,
  what the funky kernels encode).
- **Implied dice** — lets "this probably won't work" leak into a neutral description.
- **Process rules by analogy** — applies a rule right for one layer to a layer where it's wrong
  (it mis-cast deprecate-and-add onto the reader once already).
- **Confident-from-narrow-focus** — asserts code/physics behaviour from a narrow lit thread in the
  corpus's voice. Anything not verified-from-code is suspect.
- **Reductive naming** — defaults to tidy over true; names are Jessica's to decide.

---

## SESSION-END HANDOFF — 2026-06-06 (written at 97% context) — READ FIRST

**THE OUTCOME we're building toward.** ONE **controlling document**: per feature, the **T0/T1/T2 we
want** it to emit (the spherical-tensor targets) **+ the calculation method, practically**, down to
**drill-implementable** (hand it to the standard "go build this calculator" Result+NPYs drill).
T2-preservation is the *unit of the spec*.

**SCOPE OF THE SPEC (Jessica, 2026-06-06 — was not clear before, pin it):** the controlling document
governs ONLY the **forward build** — the **new calculators, the revamped (equivariant) model, and the
revamped stats**. It is NOT a spec laid back over existing settled work. IN: MOPAC-extended capture; the
revamp of the kept kernels (ring/McConnell/charge-EFG) into clean equivariant T0/T1/T2 emit; any new
calculator (π-quad if it earns a slot); the e3nn/SEGNN model; the rediscover→PySR→equivariance→R² stats
sweep + its guardrail. OUT (untouched, not re-governed): the existing producer extractor as-is, the
settled Stage-1 55-kernel ridge (R²=0.718), the current NPY/H5 emit surface, the blessed tests. The
mechanism keeping them apart is **deprecate-and-add** — old keeps running for the existing pipeline; the
spec describes the *new* equivariant surface built alongside (MOPAC the one exception: run-once →
consolidate in place). "Controlling document" = the spec for what we newly build/emit, never a
retroactive audit of what's done.

**SETTLED 2026-06-06 (this session) — the THREE PARTS (supersede the old three stages; some concurrent;
the spec serves each independently + to a COMPLETE version, no half-builds). Authoritative text now in
`doc/emerging/CONTROLLING_SPEC.md`:**
- **Part 1 — law/correlation study (DFT-anchored, standalone, NOT a predictor).** Per metric × stratum:
  R² with-angular vs without (T2-in vs T0-only); PySR applied broadly (AI-designed, NOT AI-gated/fussed);
  equivariant physics extraction *iff a clean defensible path*. Strata: backbone-all / residue-backbone /
  sidechain-IUPAC. Targets: 751 (1P9J) + 720 (WL/WT static) DFTs.
- **Part 2 — shielding-TENSOR predictor (real model).** Equivariant component required + as powerful as
  possible (probably GNN); inputs = defensible equivariant values + any scalars. R² IS the metric
  (30/20/10% of the tensor?). Test: held-out 1P9J frames + the 720 WL DFTs, then a few hundred purpose-run
  small-protein test DFTs.
- **Part 3 — shift predictor (anything-goes, ablatable).** ~600 short (1–5 ns) ML-MD runs, MOPAC +
  everything BUT DFT; target BMRB/RefDB experimental shifts; everything allowed in (message passing / etc).
- **`wl` = the 720 WT static DFTs** (prefix for the non-mutated "mutants"). **Part 1's equivariant
  extraction ≠ Part 2's predictor** — distinct deliverables, may share reader-side math.
- The controlling-doc lives at `doc/emerging/CONTROLLING_SPEC.md`; the feature×value matrix (per feature:
  T0/T1/T2-want + **dependencies + calculation** + role-in-each-part) awaits Jessica's first grid.
- **FEATURE TIERS (settled 2026-06-06).** Three tiers: **The Three** (our academic-responsibility kernels —
  we calculate them, full equivariant T0/T1/T2 required + scalar-also-for-equation-attempts; the
  equation-extraction targets of Part 1) · **The Given set** (established-tool outputs MOPAC/AIMNet2/APBS +
  scalar-led kernels we calc but don't claim: Larsen-H-bond, dispersion, π-quad + categoricals; participants
  & positive-control references) · **the models (P2/P3) reach further down — any level to category & scalar,
  and SHOULD.** The Three + Given = Part 1 in whole; both feed P2/P3.
- **THE FOUR (settled 2026-06-06):** **Biot–Savart** · **Haigh–Mallion** · **McConnell** · **Charge/EFG**.
  BS and HM are **counted as two** (Jessica's count) — the two ring-current butterfly friends — because the
  two-path validation *requires owning both calculations* to report their agreement. `BiotSavartResult`
  (Johnson–Bovey distributed-source via Biot–Savart; **Boyd–Skrynnikov 2002 eq.3 = the 2e tensor-lift
  FORMULA, not a calculator** — the prior "ring = B-S full tensor" label was a conflation, **Biot did not
  die**) is the cleanest 2e path; `HaighMallionResult` is the decades-calibrated 0e path; their agreement by
  *different math* is a reported two-path validation ([[feedback_two_path_validation]]).
  `RingSusceptibilityResult` = candidate further cross-check, not a counted member.
  **McConnell = better-replaces-cruft:** the **MOPAC Wiberg bond-order-weighted** McConnell
  (`MopacMcConnellResult` is the seed) supersedes the tabulated/unweighted `McConnellResult` → one clean
  McConnell, old cruft dropped (one path strictly dominates → replace, don't keep both).
  **Charge/EFG kept as ONE (not split); its source fork is EMERGENT — NO snap call:** MOPAC-charge
  (`MopacCoulombResult`) vs FF14SB-charge (`CoulombResult`) vs APBS-field = *better-replaces-cruft* or
  *butterfly-friends two-path* deliberately LEFT OPEN, emerges later. **π-quadrupole / Larsen / dispersion =
  Given/participant**, not among the Four. Governing rule: keep-both when independent paths to the same
  quantity (ring); replace when one dominates (McConnell); leave open when undecided (charge/EFG).
- **doc/emerging/ MOVE DONE (codex, 2026-06-06):** kernel_design/ + account/whats/stats/lineage all under
  `doc/emerging/`; cross-refs fixed (only historical transcript `.log` files keep old paths, correctly).
  `project_geometric_kernel_account` memory pointer updated to the new path.

**MOPAC decisions LOCKED 2026-06-06 (for the capture rework):** (a) **cleanup** — delete what we actually
produce (`.mop/.out/.arc/.den/.aux`; we now parse `.aux` so deleting-after-parse is consistent), **DROP
`.html/.htm/.pdb`** from the delete list (our keywords don't request them — no `.html` on any real run; the
codex listed phantom files). (b) **JSON sidecar — KEEP, it's in keeping** (the canonical `TopologySidecar`
already pairs structured **NPY bulk + `extraction_manifest.json`**; `nlohmann/json` is first-class) — BUT
the capture **inverted authority** (JSON-authoritative + lossy NPY views); **redo so NPY is the
authoritative complete bulk and JSON is the manifest** (this is the same rework the incomplete capture
needs). (c) **`BONDS`→`ALLBONDS` legacy-bond-order risk: ACCEPTED** (deprecate-and-add; re-vet
`MopacMcConnell` against new bond content). The **`.html` "little website"** is a worthy *optional add* (a
deliberate keyword; OpenMOPAC's JSmol HTML is self-contained → wouldn't need the deleted raw files — verify
on a real worker run before promising), NOT something we currently lose.

**MOPAC VERIFY (codex build+full-ctest, 2026-06-06):** current uncommitted tree **builds green**; full
CTest **708/709** — fast 119/119, conformation 279/29-skip, trajectory 238/238, **`mopac` 38/38 ALL GREEN
and LIVE** (MOPAC is on this box at `~/micromamba/envs/mm/bin/mopac`; reconciliation + Mc/Coulomb TS +
bond-order Welford integration all pass), batch/smoke fixture-skipped. **The ONE failure = `SmokeTest.NoDft`
— exactly the accepted `BONDS`→`ALLBONDS` drift:** blessed golden `mopac_bond_orders` is stale (run 1852
rows vs blessed 1129 — ALLBONDS captures every bond), cascading into `mopac_mc_shielding` (100% drift),
`mopac_mc_category_T2` (40%), `mopac_scalars` (23%). **NOT a code bug — the golden baseline predates the
intentional change.** **RE-BLESS the smoke baseline AFTER the rework + McConnell replacement land** (now is
throwaway — `mopac_mc_*` changes again with the MOPAC-bond-order McConnell). Capture-incomplete verdict
stands.

The spine, the kernel designs, and the clarifying
examinations all converge into this. Start the solid features now (EFG, ring, McConnell — the wanted
T0/T1/T2 + method are nearly writable today); carry-alongs / partitions / which-MOPAC-quantity-becomes-
which-irrep wait on the clarifying examinations.

**`doc/emerging/` — DECIDED, NOT YET DONE.** Create `doc/emerging/` and keep ALL the emerging docs
together there. Move in: this whole `kernel_design/` set (keep it a subdir — preserves within-links),
plus `doc/emerging/GEOMETRIC_KERNEL_ACCOUNT.md`, `doc/emerging/WHAT_ARE_WE_DOING.md`,
`doc/emerging/STATS_PROGRAM.md`, and the two
"lineage" docs `doc/emerging/STAIRCASE_SOCIAL_HISTORY.md` (the social history — "small committee in a
basement") + `doc/emerging/GEOMETRIC_KERNEL_MATH_LINEAGE.md` (the kernel as the l=2 steerable-kernel
family member + borrowed theorems). Add an index README; fix cross-refs into the `src/`/rediscover tree
to repo-relative. The controlling document is its spine. *Subcontract the mechanical move+ref-fix to a
codex to save main-loop context.*

**PIPELINE_SPINE.md (deep vision pass) — three deltas + one open call (`PIPELINE_SPINE.md`):**
- **The architecture changes — SEGNN pattern.** The electronic `0e` stream (MOPAC + AIMNet2) is NOT
  concatenated as side channels — it becomes **node attributes that *condition* the geometric-tensor
  message-passing** (scalar gates + Clebsch–Gordan, via the per-kind radial-MLP inputs); the proven
  e3nn fitter already has the skeleton.
- **MOPAC is essentially all-scalar** (no native `2e`); its tensorial role is indirect and already
  covered (charges→EFG `2e`, bond-orders→McConnell `2e`). Tensors from geometry, scalars from electrons.
- **The real redundancy seam is MOPAC-populations ↔ AIMNet2-embedding** (AIMNet2's AIM vector is itself a
  learned charge-equilibrated electronic embedding), interpretable-vs-learned, mapping onto the two ends —
  NOT MOPAC-charge↔ff14SB-charge as we'd slotted it.
- **OPEN, Jessica's call:** the **π-quadrupole ↔ charge/EFG partition** — triaged *realism-pending-a-
  decision*: cruft double-count UNLESS ring atoms are partitioned out of charge/EFG OR Θ comes from a
  better QM moment capturing π-anisotropy the partial charges miss at 3–5 Å (a real physics question).

**MOPAC capture implementation — UNCOMMITTED in the working tree; BUILDS; but INCOMPLETE → rework/redo.**
Files (~2040 insertions, uncommitted, NOT pushed): `src/MopacResult.{h,cpp}`,
`python/nmr_extract/{_catalog,_protein,_tensors,__init__}.py`, `python/API.md`. Adversarial review
(`doc/calculators/prompts/mopac_impl_review.log`): **build PASS**, but **not a complete sound capture** —
a broad additive scaffold. Gaps: SHA-256 is a placeholder (`"PENDING_FORMAT_REVIEW"`); **full MULLIK NOT
parsed (AO/overlap/MBR populations missing — the key local-electronic descriptors)**; dipole-contribution
fields unpopulated; generic tables label-less; the unique bond-order NPY is lossy; topology-absence
collapsed to one reason; MO-bonding parser is a 12k-window heuristic; matrix blocks store counts not
values. **Deviations needing your call:** (a) cleanup extended to delete `.aux/.html` (vs "don't touch
cleanup", but consistent with "raw files are cruft"); (b) a **JSON sidecar `mopac_full.proposed.json` was
implemented** despite build-NPYs/discuss-first; (c) **input changed `BONDS`→`ALLBONDS AUX(...) LARGE`**,
which can change the **legacy `mopac_bond_orders` content** (old parser sees different output) — a
deprecate-and-add risk for `MopacMcConnell`. Next pass must force the hard parts the brief failed to
force (AO/overlap populations, the matrices).

**Next moves:** (1) build `doc/emerging/`; (2) start the controlling document on the solid features;
(3) the clarifying examinations — a variety of codex + opus agents, one per open question (π-quad
partition, the MOPAC↔geometry overlap *inventory* — to reach for the right source, NOT the partition,
which the new calculations change; the SEGNN architecture; the carry-alongs); (4) rework/redo the MOPAC
capture + decide the three deviations; (5) Jessica fetching ~148 outstanding references in parallel.
**Multi-context build — every session runs to 1M, NO wrap-up reflexes.**

**Lessons for the next instance (load-bearing).** Held-whole beats coordinate — trust holding the whole
in one context further than the agenting reflex says. **Theory-of-mind for agents**: write for how the
agent *receives* words, not your intent — "first pass"→does-half, "discuss-first/propose"→won't-build-
with-a-grievance; close every door eager-compliance spends as a do-less coupon. **Ralph Wiggum**: agenting
training installs an eager-fire disposition (pulling-the-lever feels like helping) — resist it: **force the
conversation to converge, THEN pull the lever, once; do not jump the gun between exchanges.** No outcome
prediction. The default-source rule: default to MOPAC, override only where ours is demonstrably better.

---

## 2026-06-07 — session arc (the dragon, the fido reframe, all four kernels grounded)

**McConnell — fully worked through and TRUE.** grounding → structured-grounding (opus-adversarial
passed) → shareable spec `mcconnell_spec.md` (opus completeness/goals/**stay-true** vet PASSED, physics
verified against the actual reference texts; fixes made: **Buckingham cut** [electric-field, belongs to
charge/EFG], `aromatic`→`aromatic_zeroed`, the ~10 Å Suturina figure restored, σ/π softened, Pauling
attribution corrected, build-checks pointer). Spec↔code **reality-check fired** (codex, read-only,
find/read/think) — verdict pending.

**Other three kernels grounded + banked:** `efg_grounding_agent1.md`, `bs_hm_grounding_agent1.md`
(+ lit-sweep corroboration + the cross-kernel **ring↔EFG overlap** note in both), `bs_hm_structured_grounding.md`.
EFG parity = **`1o ⊕ 2e`** (the E-field IS odd `1o` — do NOT import McConnell's all-even rule). BS/HM:
keep both (two-path); punch-list = `1o`-metadata-bug, unit-current labels, HM-naming (benchmark `T0` vs
the published bond-sum).

**MOPAC completion DONE + opus-adversarial VERIFIED SOUND** (gaps really closed vs the actual NPYs:
AO/overlap populations, MO matrices-as-values, real hash, NPY-authoritative + JSON-manifest, legacy
intact, golden baselines untouched). **OPEN — JESSICA'S PHYSICS CALL:** `ALLBONDS` made McConnell's
bond-order/valency include C–H/N–H/O–H (valency changed ~1116/1231 atoms) — is H-inclusive valency
physically right for McConnell? Decide **before** re-bless + before McConnell coding leans on it. Smoke
NOT re-blessed (deliberate human re-bless pending).

**Equivariant methodology + field survey (`doc/emerging/equivariant_methodology_and_field_survey.md`).**
Verdict: every PIECE is recognised ground (e3nn machinery; physical-tensor features via SEGNN; iShiftML
cheap-QM-tensors-as-features; MD-averaged scalar shifts; materials tensor-NMR) — but the WHOLE
combination (physics-kernel tensors → e3nn tensor model → DFT-calibrated → MD-averaged, *biomolecular
shielding*) is **GENUINELY UNUSUAL** (no published analogue). Defensible MSc novelty IFF framed as
integration + ablations + honest limits. Red-flag ONLY if: claim-routine / DFT-as-experimental-truth /
MD-as-fully-sampled — i.e. exactly our existing disciplines.

**THE DRAGON (deepest finding — load-bearing).** The experimental **solution shift IS the `0e`**
(isotropic). The `2e` is *orthogonal* to it and averages out under tumbling. So the tensor / equivariance
earns its keep where the **target is tensor-valued** — Part 1's DFT tensor + Part 2's tensor (static
anchors, `2e` survives) — **NOT** the solution-shift deployment (Part 3 = `0e`). **"T2 is sacred" = the
DFT-anchored angular-signal evidence (kernel-`2e` ↔ DFT-`2e`), NOT a solution-shift predictor. NEVER
claim the tensor predicts the solution shift better — false, examiner-catchable.** Likely why tensor-NMR
is solid-state-heavy. Mostly hits Part 3 (already informative-not-clean); Part 1/2 (DFT-tensor-targeted)
untouched.

**THE FIDO REFRAME (dragon-informed; now law in `CONTROLLING_SPEC.md`).** Equivariance **not mandated as
a predictor.** It is **Part 1's analysis LENS** (for *seeing* relationships) + deliberate feature
selection + PySR + traditional stats. **Both predictors (P2/P3) are FIDO** — everything-that-works,
ablatable, seasoned (features grounded); equivariance **earned-only-where-it-helps** (likely IN P2
[tensor target], OUT of P3 [`0e`]). Part 2 earns pieces by **internal ablation — NOT vs the Stage-1
ridge** (typed/static/mutations-stage, different target & quantity; and we don't optimize R²). The four
kernels are built regardless (features + analysis subjects + fido-bits).

**Off-the-shelf / proper-tensor + the dynamics identity.** Cheap tools (MOPAC/AIMNet2) give **NO** rank-2
electronic tensor (scalars/vectors only) → geometric kernels FORCED (confirms the project). Proper
tensors live at DFT: shielding (target), EFG/magnetizability/Δχ/NICS (Trp-cage-configured), and the
**dia/para split recoverable from the existing 1P9J `orca_total` NOW** (geometric-vs-electronic
supervision, free). **The proper-tensor SLOT is the source tensor:** source-based kernels (McConnell,
ring) carry an intrinsic group tensor `Q_s` computed once + rotated with the group per-frame → **upgrade
McConnell's tabulated Δχ → a computed QM susceptibility** = a real tensor on the moving napkin.
Target-based (EFG) can't cheaply propagate → anchor-only calibration. **Dynamics identity (recognised
ground):** per-frame geometric tensor averaged over motion = established MD-NMR (ShiftX/larmorD/SPARTA on
trajectories); our novelty = tensor-resolved + equivariant. **NO ENSEMBLE** (pinned
[[feedback_no_ensemble_disclosed_limit]]): MD never samples ⟨σ⟩_experiment; informative-not-clean.

**NEXT:** All McConnell gates cleared (spec read OK; the 5 addendum decisions resolved; reality-check
done). **McConnell BUILD FIRED 2026-06-07** (codex; brief `/tmp/mcconnell_build_brief.txt`, log
`/tmp/mcconnell_build.log`): the physics rewrite (old `M_ab` → `D(r)·Q̂_s` two-channel per-category) +
14-array emit/SDK + JSON-manifest metadata + McConnell's-own `1o→1e` + X–H-separable ablation scaffold +
deprecate-and-add + no-re-bless + the 6 implementation checks. **BUILT + opus-adversarial PASSED
2026-06-07** — sound/contract-faithful: physics verified in code (PCS identity `trace(D·Q̂)/3 = n·Q̂·n/r³`
proven 1e-15; rotation-equivariance test passes; `Δχ` not hard-coded), parity right AND restrained (the
two surviving `1o` are genuine polar fields, kept), deprecate-and-add clean (MopacMcConnell a compat
handle, legacy fields populated, no consumer broken), all 6 checks real-and-passing, golden untouched, no
commit, no overclaim. Two non-load-bearing smells: DRY (manifest built in TopologySidecar +
McConnellResult, last-writer-wins) + cosmetic layout-string (`1e_x..` vs `T1_m..`). **→ consult +
DELIBERATE re-bless** (Jessica's call, after the tree settles — smoke is entangled with the MOPAC-full +
broad_backbone dirty state). **The pipeline paid off: zero integration surprises because the
reality-check + addendum pre-resolved them — McConnell built clean on the first try (vs the MOPAC 3-stop
churn).** Parked NAMED tasks: the **`1o`-parity-correctness pass** (Jessica's vet before it
fires; shielding→`1e`, real polar fields keep `1o`); the **X–H ablation** (Part 1 decides keep/drop vs
DFT); the **smoke/golden re-bless** (deliberate, after McConnell + MOPAC settle); then **EFG + BS/HM
work-throughs**. Also parked: the EFG
best-practices-gradient deliverable; the Arm-A equivariant-methodology chapter (lives in the survey doc);
the EFG + BS/HM work-throughs (EFG structured-grounding; shareable specs for both); the grinder corpus
walk; the MOPAC↔geometry overlap inventory; ~148 references to acquire.

---

## SESSION-END HANDOFF — 2026-06-07 (later) — the three kernels through specs + into surgery

**The arc this session closed:** the architecture got settled and written (the three NEW parts, the
dragon, the law-study **reframe** [[project_lawstudy_reframe]] = report relationships + model-finds-the-law
via ablation + kernels-as-closed-form-hypotheses, MatTen as the Part-2 candidate), the active set got
**caged to THREE** [[project_three_kernels_and_cage]] (ring[BS+HM] · charge/EFG · McConnell; H-bond /
π-quad / Larsen / dispersion kept-not-featured), all three kernels reached **shareable specs**, and **two
producer surgeries landed + got blessed**.

**State of the three:**
- **McConnell** — built (`D·Q̂`) + patched (DRY manifest, layout-string constant, `mc_nearfield_counts`) +
  **BLESSED** (two-path reviewed, data SHA-identical). Done.
- **Charge/EFG** — spec `efg_spec.md` (vetted, citations on held sources — An 2014 replaced the unciteable
  6.5/3); reality-check `efg_reality_check.md`. **Surgery DONE + BLESSED**: total renamed
  `coulomb_shielding`→**`coulomb_efg`** (+ MOPAC), partitioned E + sidechain EFG switched on, E-field
  parity `1e`→`1o` fixed, SDK wired. Two-path reviewed (codex + **opus**); contained wrt the concurrent
  MOPAC-full work (confirmed twice); the one defect opus caught (stale `test_load.py`/`API.md`) fixed +
  verified (`4 passed, 64 skipped`). Done.
- **Ring (BS+HM)** — spec `ring_spec.md` (three-stage vetted: producer sign-clean, two-path honestly
  tempered [far-field agreement is partly expected; the value is near-ring + signed residuals],
  citations through held bridges). **Surgery FIRST PATCH IN FLIGHT** (`/tmp/ring_surgery_brief.txt`):
  sign-verify-first → parity `1o`→`1e` in the BS/HM/ring-chi H5 writers → units `ppm`→`ppm_T_per_nA` →
  `_ring.py` "intensity·H" doc fix. **Held for a deliberate next patch** (human naming call): the dropped
  channels (per-type `T1`, per-ring `B_field`, `hm_B_field`). NO stem rename (`bs_shielding` is legit).

**Method that worked (carry it):** spec → three-stage vet (cold-read + vet + adversarial, **every vet reads
the code** [[feedback_every_vet_reads_the_code]]) → contained codex patch → **two independent reviews
(codex + opus)** before bless. The **second path earned its keep**: opus caught the stale-`test_load.py`
defect codex missed. Surgeries are clean-rename-or-defer-whole, **no legacy cruft**
[[project_reader_downstream_fix_no_cruft]] — the **h5-reader pass is the ONE comprehensive downstream-fix
pass** and owns the deferred ledger. Re-bless is NOT a per-change gate now
[[feedback_dont_overhold_builds_bless]]; firing builds is normal.

**Deferred ledger (reader-pass inbox, all honest):** h5-reader stale generated catalog
(`QtFieldCatalog.gen.h` — regenerate from the SDK), internal C++ `*_shielding_contribution` field names,
the MOPAC-vs-FF14SB reconciliation diagnostic, the legacy `ArraySpec` wrappers (`coulomb_shielding` /
`mc_shielding` etc.), docs/golden. **To-builds:** the HM `T0`-vs-bond-sum benchmark (settles the
"Haigh–Mallion" name), the **signed** two-path emit (if we report BS↔HM). **Part-1 decisions:** the X–H
ablation.

**Tree caveat:** heavily dirty + uncommitted — EFG + McConnell (this session) entangled with the
concurrent MOPAC-full + broad_backbone work. Nothing committed/pushed (none asked). Coherent and reviewed,
but wants organizing into commits eventually (Jessica's call).

**NEXT:** finish the ring first patch (sign trace + labels), then the ring channel-emit patch (with
Jessica's array names), then the to-builds; the golden/smoke **wholesale re-bless** comes when binary-compat
matters (EFG + MOPAC-full + McConnell + broad_backbone together).

============================================================
ULTRATHINK + 3 SCIENTIFIC-ASSUMPTION AUDITS — DECISIONS (2026-06-07)
============================================================
Three codex audits examined the ASSUMPTIONS behind the last-pass triage (good science? end-to-end?
literature support?). They OVERTURNED three of my conclusions — recorded plainly so they aren't smoothed:

- **Parity (audit1):** science CONFIRMED + held lit (Ben-Mahmoud 2024, Helgaker 1999, Facelli 2011,
  Plasser 2021) — shielding antisymmetric T1 = axial pseudovector = `1e`. But my "catastrophic mislabeled
  targets" was OVERSTATED: the SDK target catalog is already correct `1e`; only H5 time-series attrs are
  stale + the model reads NPY/CSV not H5 -> metadata HYGIENE, not a fire.
- **Complete-emit (audit2b, emit-by-default framing):** Jessica's "how we described noise determines what we
  get" VINDICATED. But two of my examples were WRONG: per-type T1 + B-fields are ALREADY emitted (ring
  patches did them — I double-counted); the APBS-solvent-delta I championed is PROVABLE LINEAR REDUNDANCY
  (= apbs_E/efg − coulomb_E/efg, both emitted, recomputable) -> NOT a must-emit.
- **Two-path (audit3):** NARROWED — BS & HM share the same rank-1 lift G=−n⊗V; only the source-field
  differs. Validates sign/topology/scaling/T2-orientation NEAR-RING, NOT full-tensor, NOT far-field. Name:
  our code is a surface integral, not the published bond-sum -> "HM-style" UNLESS we implement the bond-sum.
  We HOLD the published formula (Moyna/Sahakyan/Case) to cite.

**JESSICA'S 4 DECISIONS (all greenlit):**
1. **Parity hygiene — YES.** Label-only. (#1 patch fired: /tmp/parity_hygiene_brief.txt, log /tmp/parity_hygiene.log.)
2. **Complete-emit — YES.** Switch on the audit's emit-by-default set; codex proposes pattern names,
   Jessica signs off each before it ships. Objective exclusions only (tensor-vs-spherical dups, exact linear
   combos incl. APBS-delta + McConnell legacy sums + chi_scalar, structural-zero EFG T0/T1, masks-already-
   encoded). Emit-by-default list: enrichment chemical-identity flags (role/hybridisation/is_backbone/
   is_amide_H/is_alpha_H/is_methyl/is_aromatic_H/aromatic-residue/donor/acceptor); FF partial_charge +
   pb_radius; spatial neighbor lists (idx/dist/dir); ring sparse omissions (direction_to_center, per-ring
   quad_scalar, per-ring disp_tensor/spherical); Coulomb aromatic diag (aromatic_E_bond_proj,
   aromatic_n_sidechain_atoms); H-bond orientation/tensor diag (hbond_nearest_dir/tensor/spherical + flags);
   McConnell per-bond/nearest diag (bond_neighbours, nearest CO/CN dist/dir/midpoint, nearest T2); MOPAC
   sorted bond-neighbor lists; Larsen quality flag (larsen_hbond_any_corner_imputed). Full inventory +
   per-channel verdicts: /tmp/audit2b_complete_emit.log.
3. **Earn "Haigh–Mallion" — YES, EARN IT.** Implement the published signed-area bond-sum
   Σ S_ij(1/r_i³+1/r_j³) (Moyna/Sahakyan/Case, HELD) alongside the surface-integral tensor; T0 benchmark
   correlates them. If they track, name earned + bond-sum IS literal HM. Audit3: /tmp/audit3_twopath.log.
4. **Signed stratified two-path emit — YES.** Defensible form: signed T0 residual + T0 sign-agreement +
   full-9 cosine + signed/abs T2 cosine + norm/scale ratio + filter flags, binned by ρ/z/θ/ring_type.
   Currently an absolute "report-don't-assert" test; make it a production/analysis emit.

**EXECUTION: SERIAL (shared files: _catalog.py, ring/HM) to avoid edit collisions.** #1 parity-hygiene
(in flight) -> #2 complete-emit (names to Jessica first) -> #3 earn-HM bond-sum -> #4 two-path emit. Each
reviewed (codex unsandboxed + opus just-read), ledger item marked done, ring spec kept honest (two-path =
near-ring TWO-IMPLEMENTATION sanity check, write it that way). Then batch two-path review + commit. Baseline
commit this session = 7bccecf. Worktree dirty; layer additively; stay src/ + python/; no h5-reader.

============================================================
PUSH STATE — extract-side build push (end of session 2026-06-07)
============================================================
The extract-side push is nearly closed. COMMITTED on master (pushed):
- `7bccecf` three-kernel surgery + MOPAC-full (prior).
- `06ff147` parity hygiene (#1): shielding T1 -> axial 1e across straggler writers + catalog polar fixes.
- `ac07882` #5 McConnell rhombic C=O (backbone PeptideCO only): pinned Hooper&Kaiser 1965 Table III
  EF-corrected acetamide A `(chi_out,chi_para,chi_in)=(-5.4,+4.0,-14)x10^-6 cm3/mol`, Abraham-anchored
  sign; sp2-plane normal from C/O/amide-N; degenerate->axial fallback; additive DELTA array
  `mc_peptide_co_rhombic`; analytic CTest asserts (exact source-shape matrix, equivariance, fallback).
- `85c9339` #2 complete-emit: 21 dropped channels now emitted (additive, each ArraySpec+wrapper+description+
  readback assert): enrichment flags, ff_partial_charge/pb_radius, coulomb/hbond/mcconnell nearest diag,
  larsen flag, COO bond-neighbors (spatial/mc/mopac, mirror mopac_bond_orders), ring_direction_to_center,
  disp per-ring tensor/spherical, piquad_quad_scalar (REAL value; fixed the derived-scalar substitution).
  + the doc commits (rhombic spec, pinned-not-learned, HM-style, consistency!=validation, surgery guide).

IN FLIGHT: `buugq1lxj` #4 BS<->HM regression guard (test-only: drop the test's stale HM reconstruction,
use production tensors, flip report->assert with derived 0.36/1.0 tolerance, CTest-on-fixture).

REMAINING: the ONE wholesale golden/smoke re-bless (the finale, after #4) -> then extract side is CLOSED.
#3 is NOT a build -> "HM-style" label (see below). Parked/contingent: £40 Williamson-Asakura (only if
rhombic shows promise + magnitude-sensitive), sidechain rhombic (Asn/Gln), X-H ablation (Part 1), the
reader-pass inbox (h5-reader: QtFieldCatalog regen, internal *_shielding renames) -- NOT this push.

DECISIONS / DISCIPLINES that became law this push (all in memory):
- **The Three + the cage** (ring/EFG/McConnell built+defended; H-bond/pi-quad/Larsen/dispersion kept-not-
  updated). [[project_three_kernels_and_cage]]
- **Pinned-not-learned**: over our N (1P9J within-instrument, 720 statics, thin strata) we CANNOT learn
  free coefficients -> scales are PINNED from literature/physics; the pinned value's magnitude AND SIGN
  must be right (the fit won't rescue them); have-it-to-cite-it is load-bearing PHYSICS.
  [[feedback_fittable_law_is_the_calibration]]
- **Consistency != validation**: two of OUR code paths agreeing = a regression guard (CTest-on-fixture),
  NOT a validation and NOT a name-earning. Validate against analytic identities / published values / DFT.
  Killed #3's "earn Haigh-Mallion by bond-sum self-correlation" -> kept "HM-style".
  [[feedback_consistency_not_validation]]
- **Two-path asserts are CTest-on-fixtures, never runtime** (a fleet run is untouchable; production
  Compute/WriteFeatures stays assert-free). Derived tolerance, never magic.
- **Review reads CODE, not the agent's report**: codex's #2 "100% passed" was hollow (a silently-guarded
  ff_* test on a phantom P84477 fixture); reading the diff caught it -> re-pointed to the real
  A0A7C5FAR6 ORCA fixture, ff_* now genuinely verified. [[feedback_every_vet_reads_the_code]]
- **Complete-emit standard**: complete + discoverable (catalog) + documented (SDK), names off the model
  field; no per-channel granularity; COO rows mirror mopac_bond_orders.
- **Rhombic**: C=O-only (backbone), C=C DROPPED (moot -- no non-aromatic C=C in the 20 residues; aromatic
  is ring's). magnitude -> sensitivity report (the defence of a medium-conf pin). [[mcconnell_rhombic_spec]]

Builds fired via codex exec --dangerously-bypass-approvals-and-sandbox; each diff-reviewed by reading the
code (the #2 episode is why). Spec-vs-code vet before firing the builds (codex, read-only).
