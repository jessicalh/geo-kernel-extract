# Fix to vision — rediscover, 2026-06-02

A clear-eyed understanding pass after ~2.5h of work drifted off the founding
vision. READ-ONLY synthesis. No code changed; no prescription rubber-stamped.
Every claim is grounded in a doc or `file:line`. Where the artifacts are silent
or the vision is genuinely open, this says so. The IUPAC/topology episode is a
known trap and is not touched.

The honest headline up front: **this is not a "move a few reducers" cleanup.**
The cardinal sin (Python recompute / second model) is contained — that is real
and worth saying. But the engine the vision rests on (`one canonical traversal`)
was only *half* built by #29, the substrate the vision calls *lean and
method-agnostic* grew into a per-source × all-atom × all-frame raw dump that
Python then re-reduces and re-selects over, and the whole CSV/NPV dump sits on
the wrong side of the flagged reader-as-platform fork. Returning to vision is a
real re-shaping of the C++ emit, not a patch. The good news is that the typed
spine and the e3nn discipline are intact, so the re-shaping is *additive on a
sound foundation*, not a rewrite.

---

## A. THE ORIGINAL VISION, reconstructed from the artifacts

### A1. The spine — exactly one typed model, the C++ side computes-and-emits

`analysis/PATTERNS.md:11-19` is the constitution and it is categorical:

> "The model is the spine — ONE typed model, no second model in Python. There is
> exactly one model of the protein in this system: the typed, resident object
> model with fast lookups (`QtProtein` / topology / the `QtAtom` & `QtRing`
> virtuals + the indexed trajectory + the catalog) ... The C++ producer computes
> AND EMITS the physics off that one model — kernels, fields, sums, the spherical
> (T0/T1/T2) decomposition, local frames, the DFT target. Python never recomputes
> any of it: a recompute ... is a SECOND, shadow model of the protein living
> beside the one we define — the thing we refuse, not a style preference."

`GUIDANCE.md:22-34` makes the same point operationally and — crucially — names
the failure mode the drift would later walk into:

> "when a new analysis or relationship needs a quantity the substrate does not
> yet carry, the fix is to **extend the C++ emit** (the spine), not to compute it
> in Python."

`h5-reader/CLAUDE.md` ("The model is the spine") confirms this is project law,
not a local rediscover convention.

### A2. Emit is not a limiter — emit the LITERATURE-SCALED PREDICTION, not scratch

`EMIT_SURFACE_AUDIT.md:12-19` states the ideal precisely, with the ring as the
worked exemplar:

> "the spine asks IT the physics question and emits the literature-scaled
> PREDICTION; Python only correlates. RING is the exemplar done right — per-ring-
> type `LiteratureIntensity()` on the typed `QtRing` → emitted `jb_T*`
> (literature-scaled ppm) → `ring_literature_decirc.py` only correlates it vs
> DFT. The counter-example is the per-category geometric scalar sum that a Python
> script weights into physics."

So the vision's emit unit is `jb_T0` / `jb_T2_local_*` (RingCurrentNeighborhood
.cpp:271-272) — a typed object answering "what is your literature-scaled
shielding contribution," NOT `sum_dipolar_ringtype_*` (raw `Σ(3cos²θ−1)/r³` per
type, RingCurrentNeighborhood.cpp:280-283) that exists only for a Python OLS to
weight. `REDISCOVERY_MAP.md` step 2 names this "per-cell emit-extension … the
spine (C++) computes the physics and emits it … NEVER a Python recompute," with
McConnell explicitly required to *mirror* the ring (`#34` literature-kernel T2).

### A3. The goal — recover classical kernels from DFT, two depths, correlate-not-match

`REDISCOVERY_MAP.md:7-24`:

> "The goal: recover the classical NMR shielding kernels from DFT … Not
> prediction, not R²-optimisation: evidence the kernel set is complete enough."
> Depth A — "a law falls out (symbolic) … validated NON-circularly … the equation
> came out, not imposed." Depth B — "the signal is captured (correlate) … a
> credible per-stratum R²."

And the law/model distinction the memory `feedback_fittable_law_is_the_calibration`
encodes: the **fittable LAW** (de-circularisation, γ→1) is the only ACTUAL
calibration (arc layer 2); the **ensemble** is a MODEL (signal capture, R², arc
layer 3), not a calibration. STATE.md's status grid (line 27-38) tracks each cell
against both depths.

### A4. The engine — one canonical traversal, carrier the only thing that varies

`analysis/PATTERNS.md:20-27`:

> "The C++ functional interface is the producer-side face of that same boundary
> … express each cell through it, extend it cleanly, **don't bypass it or let it
> sprawl into one-off walks or sibling runners (the broad case forced one — #29
> unifies them)**. A blurred functional interface is a blurred model boundary —
> that is how a second model creeps in."

The engine's own header states the target shape (`RelationshipEngine.h:34-53`):
the traversal "is the ONE thing that does not vary between drivers; the CARRIER
is the only thing that does … THE REDUCER IS A PROPERTY OF THE CARRIER … not of
the engine." Layers 1-3 are real: verbs (Verbs.h), curried Layer-2 closures
(Relationship.h), the Layer-3 loop (RelationshipEngine.cpp:10-80). Minimal-and-
clarifying, no factories/ABCs (`feedback_no_abstractions`; Relationship.h:20-24).

### A5. The substrate — LEAN, method-agnostic; equivariant = real e3nn + frozen get_C()

`DESIGN.md:28-42` and `GUIDANCE.md:55-64`: the substrate is "method-agnostic
because the fitter is open" (ridge / scalar SR / equivariant SR / e3nn sum-
pooling). The *irreducible* core is the un-summed per-source set: displacement
vectors + source identity + the DFT tensor target, in raw lab frame, because
that is what an equivariant sum-pooling fitter needs and what ridge/scalar SR can
be *reduced from*. `analysis/PATTERNS.md:32-38`: equivariant = real e3nn (never
hand-rolled SH/Wigner-D), one frozen library↔e3nn change-of-basis (`get_C()`,
`change_of_basis.py`). Note DESIGN.md calls it a "**method-agnostic substrate**"
and GUIDANCE.md "a per-(atom, frame) feature/target substrate" — *lean enough to
serve every fitter*, not *a dump of every quantity*.

### A6. The flagged END-STATE / fork — the chewer rides the LIVE model; output is secondary

`spec/reader_as_platform_2026-05-29.md:16-27` (the settled architecture, lensed
2026-05-31):

> "The chewer is **scripting against the live C++ protein model** … The chewer's
> pybind11 surface exposes that model to Python; transforms are short scripts that
> query it in its native NMR-chemistry vocabulary … and emit derived rows. Output
> format (Parquet, Postgres, Arrow IPC, static snapshots) is downstream of the
> model … Spatial structure lives in the running model; the durable output is
> whatever subset of 'answers from the model' needs to be persisted."

`h5-reader/CLAUDE.md`: "The load-bearing chewer engineering is the pybind11
binding surface, not the output sink." This is the fork the brief asks to hold
honestly: in the end-state vision, **the giant CSV/NPY dump is itself the
anti-pattern** — an auto-exporter middleman between two processes, where the
intended shape is one process holding the live model and Python querying it in
place. See §C5 for the position.

---

## B. THE DRIFT MAP — where, how, how far

Magnitude framing first: the cardinal sin is contained. The python-centric audit
(`/tmp/codex-pythoncentric-audit-out.md`, verdict) is explicit: **no** Python
opening `trajectory.h5`, **no** DFT projection recompute in model paths, **no**
hand-rolled equivariance; e3nn paths and frozen `get_C()` are in the right shape.
That is the floor holding. The drift is "**reducers, selections, and scalar-law
assembly living in Python**, over a giant raw producer dump." Below, the specific
drifts, with the irreducible-vs-drift line drawn explicitly per §A5.

### B1. The engine was only HALF unified by #29 — FOUR runners still bypass the traversal

This is the biggest under-stated drift, and the brief's central #29 question.
`ALL_ATOM_EMIT_REVIEW.md:79-83` named "three drivers"; the reality after #29 is
worse. #29 (`2675565`) genuinely folded **two** carriers onto `RunTraversal`:
`RunRelationship` (ring/mc/charge_dipole, RelationshipEngine.cpp:82-96) and
`RunBroadBackbone` (BroadBackbone.cpp:524, "the SAME … traversal as …
RunTraversal"). But it left **four** hand-written `for row : dftRows() × for atom
: atomCount()` walks outside the engine:

- `RunAllAtomEquivariantEmit` — AllAtomEquivariant.cpp:399-401 (the 68 GB
  substrate's producer). The commit message admits this verbatim: *"the all-atom
  rich-record carrier keeps its own thin walk; its record type is genuinely
  different."*
- `RunEfgPerAtomFeature` — EfgFeature.cpp:170-172.
- `RunBuckinghamEfieldPerAtomFeature` — BuckinghamEfield.cpp:140 (+ its own walk).
- `RunAimnet2PerAtomFeature` — Aimnet2Feature.cpp:188-190.

So `analysis/PATTERNS.md:20-27`'s "don't … sprawl into … sibling runners" is
violated by four siblings, not one, and #29 healed the two that were *already
reducer-carriers* while leaving the four *raw / per-atom carriers* — which are
exactly the ones that drove the fatness and the Python-reduction pressure.
**#29 is partial** (`ALL_ATOM_EMIT_REVIEW.md` verdict; STATE.md:31-33 booked it).
The all-atom bypass is not an accident of "different record type" — it is the
symptom that the carrier seam was generalized for *reducer-shaped* carriers but
not for the *raw-source* and *per-target-feature* shapes. See §C2.

### B2. Reducers living in Python — the C++ reducer done in numpy

The audit's "Python Reassembly" section, all confirmed against the scripts:

- `mcconnell_dchi_calibration.py:207-268` source-sums emitted McConnell T2 rows
  by `row_id`/`bond_category` with `np.add.at()` (`:255-256`), then deweights the
  C++ scalar (`:271-282`). This is a C++ reducer (per-category McConnell-Δχ T2)
  done in Python.
- `ensemble_layer3.py:435-470` repeats the same per-category source sum;
  `:590-593` compares the Python sum to the C++ aggregate.
- `johnson_bovey_region_recovery.py:299-320` groups emitted source `dipolar`/
  `jb_T0`/`jb_T2_local_*` by `gid`.
- `ring_literature_decirc.py:421-439` groups `jb_T*`/`jb_unit_T*` by `gid` (less
  bad — sums blessed source kernels — but still a reducer Python carries).
- `look02_self_vs_throughspace.py:26-52` ranks nearest sources, splits
  `sum_all/sum_near/sum_rest` in Python.

These exist because the spine emits the *un-summed source rows* and a *raw scalar
sum* (`sum_dipolar_*`) but NOT the C++ literature-scaled aggregate the reducer-
comment is waiting for: `look03_coefficient.py:86-87` literally says *"intensity-
weighted aggregate intentionally omitted (it is a C++ reducer to add, not a
Python re-sum of `ring_intensity*dipolar`)."* The codebase KNOWS the drift and
named the fix; the fix was not built. `EMIT_SURFACE_AUDIT.md` bucket C (C1-C7)
is this backlog enumerated column-by-column.

### B3. Selections / model-queries living in Python

The audit's "Selection As Physics," confirmed:

- Stratum maps re-derived in Python across at least seven scripts:
  `equiv_t2_backbone_e3nn.py:96-110,153-175,547-550`;
  `static_environment_calibration.py:22-29,46-60`;
  `variance_decomposition.py:22-29,46-60`; `aimnet2_ceiling_ensemble.py:26-49`;
  `ensemble_layer3.py:73-87`; `mcconnell_literature_decirc.py:358-372`. The
  stratum is a typed-model property (`QtAtom` predicates) re-implemented as a
  Python name/`frame_variant` map.
- Source-identity selection by geometric proxy: `sumpool_t0.py:37-40` uses
  `SELF_R = 3.0` as a self/bonded proxy; `look02:26-30` infers self-ring by
  nearest distance — both are model questions (the typed `is_self_or_bonded`
  flag already exists on the source rows, RecordSink.cpp:141-145) answered by a
  Python heuristic instead.
- Atom-role inferred from `atom_name`: `differencing_system_id.py:288-307`,
  `mcconnell_dchi_calibration.py:142-156`.

The audit's own framing is exactly right: this is a "role/selection model in
Python," not a full second protein model — a partial second model, which
`analysis/PATTERNS.md:11-19` still refuses on principle (a blurred boundary is how
the full one creeps in).

### B4. The substrate fatness — per-source × all-atom × all-frame raw dump

From the EMIT CODE (the 68 GB was deleted; this is reasoned from emit, not `du`):

- `AllAtomEquivariant.cpp:399-401` is `for row : dftRows() (660)` × `for atom :
  atomCount() (846)` = 558,360 target rows, then per row a multiplicity of source
  rows: rings (`:460-472`), every bond category (`:475-487`, the `AllBondMidpoints`
  cloud — *every* bond, not the McConnell subset), FF14SB + AIMNet2 charge sites
  (`:490-521`), and the per-target APBS/EFG/AIMNet2/MOPAC + 256-d embedding feature
  payloads (`:524-610`). The commit `dca30b8` records `charge_mopac 30472756` rows
  and 558,360 MOPAC-shielding rows from one re-emit.
- Sidecars: target T2 / σ_iso / raw, APBS, EFG, AIMNet2, MOPAC, 256-d embedding
  (AllAtomEquivariantSink.cpp:146-159).

The audit's split is the right one and I concur after reading the emit:
- **IRREDUCIBLE for the equivariant model** (§A5, correct, NOT drift): per-source
  rows with `disp`, source identity, orientation/axis/tensor payloads, and the
  lab-frame target T2. That row-count explosion is what an e3nn sum-pooling fitter
  consumes; baking the sum in would lock out the method the problem is shaped for
  (`DESIGN.md:28-42`).
- **DRIFT (cuttable)**: pre-combined scalar scratch (`dipolar`, `inv_r3`,
  `q_over_r3`, the `sum_*` reducers, `field_z`/`field_mag` which are literally
  `.z()`/`.norm()` of an emitted vector — BroadBackbone.h:126-127); duplicated
  target/bare payloads (EMIT_SURFACE_AUDIT.md D1-D2: DFT target emitted up to 3×
  per row-kind ×2 frames; bare kernel on both row kinds + NPY); five hand-rolled
  `writeNpyF64` + five hand-rolled CSV header blocks (D3-D5).

The audit's sober magnitude estimate, which I accept: column-width cuts on the
all-atom case save single-digit to low-tens of percent; the order-of-magnitude
shrink is on the *legacy scalar surfaces* where a raw per-source dump is replaced
by a C++ aggregate (removing the source multiplicity). So **the all-atom raw
substrate is mostly irreducible-and-correct**; the fatness war is won on the
*reducer/duplicate* columns and on *not also emitting the scalar surfaces raw*.

### B5. The dump itself is on the wrong side of the flagged fork (§A6)

This is the deepest, most uncomfortable drift, and the one the other audits do
not fully name. The reader-as-platform end-state (`spec/reader_as_platform`) is a
pybind11 chewer scripting the *live* model; the rediscover engine as built is a
**headless batch exporter** (`main_extract.cpp`, CSV + NPY to disk) whose output
Python then re-reads, re-reduces, and re-selects over. Every drift in B2/B3 is a
*direct consequence of the dump architecture*: because Python receives flat files
instead of a live model, it re-derives strata, re-sums sources, and re-infers
roles that a `frame.kernel('BS', atom).T0` / `atom.stratum` / `ring.member_atoms`
binding would have answered in place. The dump is not merely fat; it is the
*shape* that manufactures the reducer-and-selection drift. Hold this honestly:
within the spike's own terms (a one-shot CSV substrate, `GUIDANCE.md:36-50`) the
dump is the agreed deliverable; measured against the *flagged end-state*, the dump
is the middleman the architecture was meant to dissolve.

### B6. When — the drift is traceable and recent

The spine/emit discipline was sound through the ring exemplar and the broad-
backbone exemplar (STATE.md:331-357, "no recompute, lead re-grep clean"). The
fatness and the four-runner sprawl entered with the *breadth* push: broad-backbone
forced the first sibling (STATE.md:439-445, "broad needed a SIBLING runner …
#29 … Sibling runner kept meanwhile"), the all-atom emit (`2ecbe59`) added the
second + the 57 GB, and the per-atom-feature cases (efg/buckingham/aimnet2) added
three more. The last ~2.5h (`2675565` #29, `dca30b8` #51, `a703701`, `cd8d710`)
unified two carriers and *welded MOPAC onto the un-unified all-atom sibling*
rather than first finishing the engine — exactly the "hard floor" fallback the
review offered (`ALL_ATOM_EMIT_REVIEW.md:143-149`), taken instead of the
preferred unification. So the drift is not ancient; it is the breadth phase
outrunning the engine, and the most recent commits chose breadth-on-the-sibling
over finishing the seam.

---

## C. WHAT "FIXING TO VISION" CONCRETELY LOOKS LIKE

Not a patch. The target shape, component by component.

### C1. What the spine EMITS — lean law-predictions + reduced kernels + model-query metadata

Three kinds of column, and *only* three:

1. **Irreducible raw per-source rows** for the equivariant fitter (KEEP, §A5):
   `disp_local/lab`, source identity (typed `QtRing`/`QtBond` virtuals),
   orientation/axis (even-form, per the `dca30b8` parity contract on
   AllAtomEquivariant.h), the DFT target tensor. This is the substrate's reason to
   exist and is correct as-is.
2. **C++ literature-scaled PREDICTIONS** (the `jb_T*` ideal, §A2): for every
   mechanism, the typed object's literature-scaled contribution AND its C++
   aggregate. Concretely the audits' convergent list:
   - aggregated ring `jb_T0` / `jb_T2_local_*` (intensity-weighted, the reducer
     `look03_coefficient.py:86-87` is waiting for) → retires `sum_dipolar_ringtype_*`
     (C1) and `sum_dipolar_all/_producer_valid` (C2).
   - per-category + summed McConnell-Δχ `mc_lit_T0`/`mc_lit_T2_local_*` (the
     machinery exists — `bondKernelT2FromSources`, BroadBackbone.cpp:79-105) →
     retires the Python `np.add.at()` reducers in `mcconnell_dchi_calibration.py`
     and `ensemble_layer3.py` (C3/C4).
   - keep raw `field_local`/`mu_local`; drop `field_z`/`field_mag` (they are
     vector accessors, C6); the Buckingham *shielding* prediction is genuinely open
     (REDISCOVERY_MAP buckingham stub ❌) — emit it when the coefficient exists,
     not a scalar scratch meanwhile.
3. **Model-query METADATA** so Python never re-selects (§B3): per-row `stratum`/
   `atom_role`, `is_self_or_bonded` (already on ring source rows — extend to all),
   source region/rank, and per-feature `support_count` / per-stratum effective-N
   in the manifest (the provenance gaps `ALL_ATOM_EMIT_REVIEW.md:240-299` and
   `cd8d710` began; finish them). Effective-N in the manifest is non-negotiable —
   computing it in Python is a recompute (`analysis/PATTERNS.md` rule 5).

Plus the pure DRY/duplicate consolidation (EMIT_SURFACE_AUDIT.md D1-D5): target-
once + bare-once (back-port the broad sink's fix to RecordSink), one `WriteNpy`,
one shared header builder.

### C2. What the ENGINE looks like — and does the sibling runner heal, or is it deeper?

The brief's sharp question: is the all-atom sibling a symptom of something
deeper? **Yes — and #29 as landed proves the seam is not yet general enough.**

#29's `RunTraversal` + `PerRecordSink` (RelationshipEngine.h:51-60) is the right
abstraction *for reducer-carriers*: it hands the carrier the fully-attached
`NeighborhoodRecord` and lets the carrier own the reduce+sink. `RunRelationship`
and `RunBroadBackbone` fit it. But it did NOT absorb:
- the **raw-source carrier** (all-atom, reducer=none) — the commit excused it as
  "genuinely different record type" (`AllAtomEquivariantSourceRecord` vs
  `SourceSlot`). That excuse is the tell: the all-atom record is richer (APBS/EFG/
  AIMNet2/MOPAC/embedding payloads), so it could not ride `NeighborhoodRecord`.
- the **per-target-feature carriers** (efg/buckingham/aimnet2) — these have *no
  source set at all*; they are pure per-(atom,frame) feature reads, which the
  source-centric traversal does not model.

So the deeper issue: `RunTraversal` is built around `NeighborhoodRecord` (a fixed
source-slot shape). The vision's "one traversal, carrier varies" needs the
*record type* to be a carrier property too — the traversal owns "for each
(atom,frame): build the local frame, the target, the source list" and hands the
carrier a *generic attached context* it interprets. The fix that *fully* heals it:

- Make the all-atom case a relationship with `reducer = none` and a raw-source
  carrier, AND generalize the source-slot to carry the rich per-source payloads
  (or let the attachers write into a carrier-owned record) — `ALL_ATOM_EMIT_REVIEW
  .md:124-134` prescription 2 is correct in spirit but underestimates the record-
  type generalization it requires.
- Model the per-target-feature cases as relationships with an **empty selector
  list** (no sources) + target-feature attachers, so efg/buckingham/aimnet2
  collapse into the one traversal too. This is what dissolves four siblings, not
  one.
- Fold the remaining triplicated geometry onto verbs — #29 did the disp/heavyParent
  fold (Verbs.h:126-149, the commit's second half); that part is genuinely done
  and good.

Net: the sibling runner is a symptom that the seam was generalized over the
*reducer* axis but not the *record-shape* and *no-source* axes. Healing #29
fully means one `RunTraversal` whose carrier owns (record type, reducer, sink),
with all six cases as carriers. That is more than #29 did, and it should be
designed across all six shapes at once (`feedback_functional_api_minimal_
clarifying_abstraction`) — not whack-a-mole, and explicitly NOT a plugin/ABC
framework (the carrier stays a `std::function` bound in code).

### C3. What PYTHON shrinks to — a thin fitter

Per `analysis/PATTERNS.md:28-31` + the audit's "Genuinely Fine": ridge/OLS, PySR,
DeepSets/e3nn training, correlations, centering, train/test splits, reporting.
Nothing else. Once C1 lands, the reducer scripts (B2) read a C++ aggregate column
instead of `np.add.at()`; the selection scripts (B3) read an emitted `stratum`/
`atom_role`/`is_self_or_bonded`/`support_count`; and the `backbone_distill_
evidence.py:257,294,329` `r^-3`/`q*r^-3` analytic comparators (the one real
recompute the audit flagged crossing the line) are deleted or replaced by reading
the emitted `dipolar`/blessed kernel. The labeled fixture (`test_change_of_basis
.py`) is the one and only surviving recompute, correctly (rule 5).

### C4. How the substrate shrinks

- Largest shrink (order of magnitude on the legacy scalar surfaces): the C++
  aggregate `jb`/`mc_lit` columns let the scalar-fit consumers stop needing raw
  per-source rows for the *reduced* mechanisms — but note the all-atom equivariant
  substrate still needs its raw rows (that is the irreducible core).
- Medium shrink on the all-atom case itself: drop the scalar scratch + duplicates
  (D1-D6, C6, C7) — single-digit to low-tens of percent of column width.
- The honest number: the substrate does not shrink to "small." It shrinks to
  "lean for what the equivariant fitter genuinely needs," which is still large
  because per-source × all-atom × all-frame is the model's actual input shape.
  The fatness *story* is "we also dumped the scalar reassembly surface and the
  duplicates," and that is what gets cut.

### C5. The fork — is the dump itself the drift? A reasoned position

**Position: for the law-discovery / equivariant work, the on-disk substrate is
NOT the drift and should stay a file dump; for the scalar/ensemble/selection work,
the dump IS the drift, and that work is what the reader-as-platform chewer was
designed to absorb.** Reasoning:

- The e3nn fitter (`equiv_t2_backbone_e3nn.py`) needs a static, shuffleable,
  epoch-iterable tensor of per-source rows + targets. A file (NPY/Parquet) is the
  right artifact for that — a live pybind11 model queried per-minibatch over many
  epochs would be slower and buys nothing, since the fit needs the *whole* tensor
  resident anyway. `reader_as_platform:23-27` agrees: "the durable output is
  whatever subset of answers needs to be persisted" — the equivariant training
  tensor is exactly such a subset. So the *equivariant* dump is a legitimate
  persisted answer, not a middleman.
- The reducer-and-selection drift (B2/B3), by contrast, is precisely the work the
  chewer dissolves: a transform that asks `atom.stratum`, `ring.member_atoms`,
  `frame.kernel('BS', atom)` of the live model and emits *already-reduced,
  already-selected* rows. Those Python scripts re-doing C++ work exist *because*
  there is no live binding — they get flat files and rebuild the model's answers.
- Therefore the dump is not monolithically right or wrong. The correct end-state
  is: the **chewer/live-model boundary** (or, near-term, fatter C++ emit columns
  per C1) answers the reduction/selection questions; the **file substrate** carries
  only the irreducible equivariant tensor + the persisted blessed predictions.
  The 68 GB symptom was the dump trying to be BOTH the equivariant tensor AND the
  raw material for every Python reduction — that double duty is the drift.

Genuinely open (do not invent): whether to build the pybind11 chewer *now* vs
keep extending the C++ emit columns as the near-term proxy. `reader_as_platform`
is "settled architecture" but the rediscover work is the never-merged spike
(`GUIDANCE.md:43-50`); the spike's contract is the CSV substrate. The principled
read is that the C1 emit-extensions ARE the chewer's `jb_T*`-style answers,
delivered as columns instead of a binding — same physics-on-the-spine boundary,
lighter mechanism. Building the full chewer is a larger commitment the lead should
decide explicitly; it is not required to *return to vision*, only to *reach the
end-state*.

---

## D. SALVAGE ASSESSMENT of the last ~2.5h

For each, KEEP / REVERT / REDO-TO-VISION with reasoning. I am willing to say
redo; nothing is preserved merely because it exists.

### D1. #29 engine-seam unification (`2675565`) — KEEP, but REDO-TO-VISION the unfinished half

The half that landed is correct and on-vision: `RunTraversal` + `PerRecordSink`,
two carriers folded, disp/heavyParent onto verbs (Verbs.h:126-149). Byte-identical
gated. KEEP that. But it is *partial* (§B1, §C2): four runners still bypass it, and
the all-atom bypass is structural, not incidental. REDO-TO-VISION the rest:
generalize the carrier to own (record-type, reducer, sink) and fold all-atom +
the three per-atom-feature cases into the one traversal. Do NOT treat #29 as done.

### D2. MOPAC-add (`dca30b8`, #51) — KEEP the wiring, ACCEPT the debt it took on

The MOPAC physics wiring is correct and on-vision: it captured the moderate
Stage-1 field/EFG leg that was on the floor (`mopacCoulombShielding()`/
`mopacMcShielding()` as T2, honestly labeled "EFG-DERIVED shielding T2, NOT raw
EFG"), with the Welford-mean static charge source and the reconciliation QC
column. The discipline held: no fourth driver (the explicit hard floor,
`ALL_ATOM_EMIT_REVIEW.md:143-149`), typed Catalog ArrayIds, provenance tags. KEEP
the leg. But it was welded onto the *un-unified* all-atom sibling — it took the
"hard floor" fallback instead of the preferred "unify first." That was a
defensible capacity call, but it means MOPAC now rides a runner that §C2 says must
be folded; when #29 is finished, the MOPAC selectors/attachers move with it. So:
KEEP the physics, BOOK that it sits on the runner slated for folding. Not a revert.

### D3. The all-atom emit (`2ecbe59`, foundation) — KEEP the raw core, REDO the fat + the bypass

Model discipline is clean (typed spine, no string-dispatch, genuinely raw lab
frame — `ALL_ATOM_EMIT_REVIEW.md:19-61`). The raw per-source equivariant rows are
the irreducible core (§A5) — KEEP. But the case (a) bypasses the engine (REDO via
§C2), (b) emits the scalar scratch + duplicate surface (REDO via C1/D-buckets),
(c) carries the provenance gaps (partly addressed by `cd8d710`). So: KEEP the raw
tensor it produces; REDO-TO-VISION its driver shape and its column surface. Do not
revert the emit — re-deriving the raw substrate is expensive and correct.

### D4. The three outstanding review bugs

- **Welford presence (`a703701`)** — KEEP. Genuine correctness fix: presence now
  gates on `n_frames_per_atom[atom] > 0 && isfinite(mean)` and `value()` defers to
  `present()` (Catalog.cpp diff). A static charge source with zero frames or NaN
  mean was previously reported present. Correct; on-vision; ship.
- **Provenance population (`cd8d710`)** — KEEP as far as it goes; INCOMPLETE.
  It populated `rowAlignmentContract` / `normalization` / `featureSupport`
  (OutputManifest.h:27-35) into the broad sink + main_extract. But the
  `ALL_ATOM_EMIT_REVIEW.md:240-299` gaps (atom-role/stratum echoed INTO the NPY
  sidecar identity, per-stratum effective-N) are only partially closed. KEEP, and
  finish the per-stratum support metadata as part of C1's model-query columns.
- **bond-axis-1o consumer** — this is a LIVE RESIDUAL, correctly handled as a
  CONTRACT not a fix (STATE.md:185-192; `dca30b8` parity §C6 on
  AllAtomEquivariant.h): raw `bond_axis_local_*` is index-oriented (producer
  min/max), sign-arbitrary, so it must be ingested as an EVEN form (outer product /
  dipolar), never as a polar 1o. The emit documents this; the e3nn consumer must
  honor it. KEEP the contract; the open item is *verifying the consumer obeys it*
  — a Python-side check, not a C++ change. Flag for the fitter author.

---

## E. THE PATH BACK — phased, gated, disk-aware

Constraints honored: no full re-emit until the emit is lean (the lead just had a
disk blowup from careless re-emits; STATE.md:48 ~30 GB free); gates are explicit
oracle parity (`analysis/oracle_parity.py --bin <h5reader_extract> --run <single .LGS>
--work <tmp> --case all`, currently ring+mc parity only) + scoped
`QT_QPA_PLATFORM=offscreen ctest --test-dir build/linux-gcc -R h5reader_rediscover`;
the extractor/library is untouchable; nmr_extract outputs under /shared are sacred.
Smallest set, in order:

**Phase 0 — finish the cheap, no-re-emit fixes (the follow-ups already started).**
Complete the provenance/effective-N manifest metadata (D4) and the bond-axis
consumer verification. These are additive, gated by ctest, and need no re-emit.
This stabilizes the foundation before reshaping.

**Phase 1 — emit the C++ blessed aggregates (C1.2), additively, alongside the
scratch.** Add aggregated ring `jb_T0`/`jb_T2_local_*` and per-category/summed
`mc_lit_*` as NEW columns; do NOT cut `sum_dipolar_*` yet. Gate: oracle parity
(the new columns are additive, ring/mc byte-identical on existing columns) +
ctest. Re-emit is required here but it is the *productive* re-emit — and it must
be drop-old, single, on the dense LGS (the `dca30b8` discipline: "old emit dropped
first, not doubled"). This is the one re-emit Phase 1-2 share; sequence the column
work so it happens once.

**Phase 2 — migrate the Python reducers/selectors to the new columns, then cut the
scratch.** Repoint the ~10 scalar-fit scripts (C2 consumer list) to the aggregated
`jb` column; repoint `mcconnell_dchi_calibration.py`/`ensemble_layer3.py` to the
`mc_lit_*` aggregates; repoint stratum/role/self-bonded selections to the emitted
metadata columns; delete the `backbone_distill_evidence.py` `r^-3` comparators.
THEN retire `sum_dipolar_*`, `field_z`/`field_mag`, and the duplicates (D1-D7).
Gate: each script reproduces its prior number off the new column (a Python-only
diff, no re-emit); the cuts are CSV-only so no SDK-contract change
(EMIT_SURFACE_AUDIT.md:405). Disk: this is where the substrate actually shrinks.

**Phase 3 — finish #29 (the engine).** Generalize the carrier to own
(record-type, reducer, sink); fold all-atom (reducer=none, rich record) and the
three per-atom-feature cases (empty selector list + target attachers) onto
`RunTraversal`; the MOPAC selectors/attachers move with the all-atom fold. Design
across all six shapes at once. Gate: oracle parity byte-identical on every case +
ctest. This is pure refactor (no physics change), so it can be gated to byte
identity and needs only a verification re-emit (drop-old), not new data.

**Phase 4 (decision, not code) — the fork.** Take the lead's explicit decision on
pybind11 chewer vs continued column-emit (§C5). Until decided, the column-emit
proxy (Phases 1-3) IS the on-vision behavior; the chewer is the end-state, not a
prerequisite for returning to vision.

Order rationale: 0 before 1 (stabilize before reshape); 1 before 2 (never cut a
column before its replacement is consumed — EMIT_SURFACE_AUDIT.md:373); 2 before 3
(shrink the substrate before refactoring the engine that emits it, so the
verification re-emit in 3 is cheap); 3 absorbs MOPAC's runner naturally. Each
phase is independently gated and (1 excepted) re-emit-free or drop-old.

---

## Bottom line

The vision is intact and recoverable: one typed spine, C++ emits literature-scaled
predictions, Python only fits, one canonical traversal, a lean method-agnostic
substrate for an open (e3nn-included) fitter, with the law-vs-model arc and the
reader-as-platform end-state behind it. The drift is real and more than cosmetic —
**#29 is half-built (four runners still bypass the engine), the Python side carries
reducers/selections that belong on the spine, and the substrate became a double-
duty raw dump that manufactured that Python drift** — but the cardinal sin is
contained and the raw equivariant core is correct. Returning to vision is a
phased, additive reshaping of the C++ emit (blessed aggregates + model-query
metadata → migrate Python → cut scratch → finish the engine), not a rewrite, and
not a quick reducer-move. The last 2.5h is mostly KEEP-the-physics / REDO-the-
shape: the MOPAC leg and the Welford fix are sound, #29's landed half is sound, but
#29 is not done and the all-atom case needs its driver and column surface redone to
vision. The fork is genuinely open and is the lead's call; it is not a blocker for
returning to vision.
