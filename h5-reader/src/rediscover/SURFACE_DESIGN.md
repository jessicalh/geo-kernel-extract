# Multi-scenario surface — design sketch (2026-05-31)

Status: **design sketch**, Step 1 of `PLAN.md` ("pseudocode the end-state").
Not built. The next session starts here rather than re-deriving. Reads with
`PLAN.md` / `STATE.md` / `FINDINGS.md` and `DESIGN.md` (the one-off). Living —
expect it to shift as the build and the codex IPC research land.

## What this is

One resident, **immutable** C++ body, loaded once; every scenario traverses it.
The one-off ring/McConnell extractions are replaced by a generic surface: a few
**primitive access verbs** on the body + per-relationship **combinator lambdas**
+ one **pure engine loop**. Python sits on top (interface + fitter-consumer);
it never re-reads the raw H5/NPY.

Shape: functional **style** (pure transforms over an immutable body, fold-based
accumulation, order-free) on a data-oriented **substrate** (contiguous arrays,
views, index arithmetic), **lazy/streaming** (no materialized cross-products, no
copies — that is what "drowning in overlapping reads" means, and we don't).

## On "pluggable interfaces" (clarified — Jessica, 2026-05-31)

"No pluggable interfaces" forbids little library **kits people trade around with
fancy registration/plugin machinery** — it does NOT mean ABC polymorphism.
- **Used / fine:** virtual dispatch on typed domain objects
  (`QtRing::LiteratureIntensity`, `QtAtom` predicates) — objects answer questions
  about themselves.
- **Fine:** the relationship combinators below — direct, **named** lambdas
  composed in code.
- **Forbidden:** a scenario plugin/registry framework, config-string factories,
  swappable third-party modules. Relationships are **named bundles defined in
  code**, not discovered or loaded.

## The supported set (agreed — Jessica, 2026-05-31)

Nine relationship TYPES; `charge_source` is a parameter axis, not separate items.
The DFT tensor (raw 3×3 + library T0/T1/T2 + local-frame) is the shared target;
each names its bare-kernel cross-check where one exists.

1. **ring_current** — aromatic C–H ← rings (`slots`, ring_nbhd) — T0+T2 [oracle; xcheck bs]
2. **mcconnell** — backbone HN ← aniso-bonds (`near`, KD) — T0+T2 [scalar done; xcheck mc]
3. **buckingham_efield** — atoms ← APBS E-field (l=1) — T0/T1
4. **efg** — atoms ← APBS EFG (T2) — T2
5. **charge_dipole** — atoms ← local Σqᵢ(rᵢ−r); `charge_source∈{ff14sb,aimnet2,mopac}` — T0/T1
6. **charge_quadrupole** — atoms ← local traceless Σqᵢ(3rr−r²I); same `charge_source` — T2
7. **larsen_hbond** — exchangeable H ← donor/acceptor geometry — T0+T2
8. **charge_response_gradient (AIMNet2 CRG)** — atoms ← AIMNet2 charge-response-gradient — T0+T2.
   NOT a Buckingham polarizability — it is `d(Σqⱼ²)/dr` (parity-odd), per
   `QtAimnet2Group.h`; do not mislabel it "polarizability." (Renamed 2026-05-31.)
9. **aimnet2_embedding** — atoms ← per-atom 256-d learned embedding (per-atom-feature kind; payload to a sidecar NPY)

Picks (settled):
- **charge_source = PARAMETER, not separate items** — the kernel is identical;
  only the charge array read differs; `substrate_conventions` already makes
  `charge_source` a required enum (no default). The FF14SB/AIMNet2/MOPAC
  sensitivity study is a batch-file param **sweep**, not N definitions.
- **embeddings IN** (#9): the equivariant machinery is built anyway, so a
  needle-moving learned-rep result is a welcome couple-paragraphs. (Structurally
  a per-atom feature, not a source sum — surfaced as an issue, not deferred.)
- **polarizability KEEP** (#8).
- **NO DEFERRAL — "there is only today"**: the design is complete now; every item
  fully stubbed + a worked example; stubs are skeletons, not punts; issues found
  are findings-to-fix, never deferrals.

MOPAC + APBS + AIMNet2 all in: APBS = PB field/EFG (#3,#4); FF14SB/MOPAC/AIMNet2 =
charge sources (#5,#6); AIMNet2 response-grad = #8; AIMNet2 embedding = #9. MOPAC
charge data lands next morning; stub + example now regardless. #5/#6 are
analysis-side multipoles off resident charges — NOT re-adding the retired
production Coulomb calculator.

## T2 frames for the new items (RESOLVED — not a punt)

T2 is the thesis; the new items' tensors are NOT deferred to "later." The
equivariant fit is what unlocks them, mostly **in the lab (MD/H5) frame, today**.
(Resolves the STUB_LANGUAGE.md L1/F2 issue and supersedes the earlier
"T2 lands lab-frame / T0 is safe" shrug.)

- **Field items (efg, efield, polarizability)** — the APBS EFG is already a T2
  computed in the MD frame, and the DFT T2 is in that same frame (verified:
  ORCA-vs-H5 rotation ~1e-4°, `CheckDftFrameAlignment`). So it is a DIRECT
  tensor-vs-tensor relationship: fit DFT-T2 from the field-T2 **equivariantly**,
  which is **tumbling-safe by construction** — both tensors rotate together under
  protein tumbling, and equivariance respects the shared rotation. **No per-atom
  frame needed.** The trap: a naive component-wise correlation in the lab frame
  IS tumbling-contaminated — the equivariant fit is precisely why you don't do
  that. (Rotation-invariants T0 / |T2| are safe either way.)
- **Source-sum T2 kernels (charge_dipole, charge_quadrupole)** — the kernel T2 is
  built from source displacement vectors; **fit equivariantly in the lab frame**
  (same tumbling-safety) — that path needs no per-atom frame and is the
  recommended route. A per-atom local frame is an *alternative* (tumbling-removed)
  representation, but **that frame does NOT yet exist**: `LocalFrameBasis` has
  only HN + aromatic-H, and the conventions doc has typed HA/Cα/CO frames, NOT a
  universal rule. So a generic per-atom frame is a convention to **DEFINE + BUILD**
  (collision-safe and TYPED, not positional — the IUPAC-trap), not an existing
  resolution. (Corrected 2026-05-31 after both dry-runs caught the overclaim.)
- **ring_current, mcconnell** keep their existing typed frames (aromatic-H / HN),
  already resolved.

Net: T2 on the new items is real now — field items via the equivariant fit in
the lab frame; source-sum kernels via equivariant-in-lab-frame or the generic
bond-local frame. This is the payoff the equivariant path was built for — which
retro-justifies building it first.

## Output carrier (DECIDED — 2026-05-31)

NOT one superset CSV (that was the one-off's shape; it doesn't scale to nine
heterogeneous bundles — the embedding's 256 columns would leak into the ring CSV).
- **Per-relationship schema**: a bundle's attachers + reducer DECLARE the columns
  they emit; `SchemaFor(relationship)` builds that relationship's schema from the
  declarations — no shared superset.
- **Two relationship kinds** under the shared `(atom,frame)` outer loop + shared
  identity + shared DFT-target columns: `relationship_kind = source_sum |
  per_atom_feature`. The reducer + row shape branch on it; per-atom-feature is
  **first-class** (read the per-atom value → row), NOT a degenerate `self()`.
- **Wide / array-shaped payloads → a new set of NPYs, documented in the SDK**
  (Jessica, 2026-05-31): the substrate emits NPYs for the array data — the 256-d
  embedding, the 5-component T2 (source and target), per-source matrices — each
  with a one-`ArraySpec`-per-NPY entry in the SDK catalog
  (`python/nmr_extract/_catalog.py` discipline; read-only wrapper). CSV carries
  identity + scalars + small features; the NPYs carry the arrays, keyed by
  (atom, frame[, source]); `manifest.json` references both. No 256-column CSVs,
  no Vec3-clobbering of a 5-vector.
- **`cross_check_kernel` is OPTIONAL** — present only where a producer bare kernel
  exists (ring=bs, mc=mc). **APBS EFG is the SOURCE for item #4, not a
  cross-check** (a field gradient in V/Å², not a ppm shielding kernel).
- **ValidateScenario phase** runs BEFORE the traversal: an absent charge_source
  (or any missing required array) fails loud there with a clear rc, never inside
  the hot loop, never a silent fallback.

## Resolved decisions (2026-05-31 chat)

**#1 — FF14SB charges: ADD THE READ TO THE READER (Jessica), contract-safe.**
De-risked against the actual `topol.top`: the protein partial charges are
**inline** in the single `[ atoms ]` section (charge = column 7), NOT behind a
molecule `#include` (the `#include`s are only `amber14sb.ff/forcefield.itp`,
`posre.itp`, water, ions). So the read needs **no `#include` resolution, no
glob**:
- a reader-side `topol.top` `[ atoms ]` parser → per-protein-atom partial charge;
- attach to an **additive `QtAtom` partial-charge field** (or a charge accessor);
- map **1:1 to protein-atom order** with a `(resnr, atom-name)` cross-check vs the
  loaded model; **fail-loud on any mismatch** (no silent reorder);
- path from the `.LGS` `topology_top_abspath` (documented convention — not
  discovery). Do NOT break the reader contract: into the typed model, additive.
- `charge_source ∈ {ff14sb (this read), aimnet2 (H5, real now), mopac (per-frame,
  writing now)}`; an absent source → **ValidateScenario fails loud**, never a
  fallback.

**#3 — typed-atom lookup contract (IUPAC-trap-safe).** Lookup is SCOPED to a
chemically-defined set (a residue's atoms or a ring's members — both already
exist) and filtered by a partial TYPED selector (any subset of Element / Locant /
BranchAddress / DiastereotopicIndex / BackboneRole; unspecified = wildcard). It
returns a **SET** (size 0/1/N); multiplicity is real chemistry (an equivalence
class), **never resolved by atom index/order/position** (that is the IUPAC-revert
trap). Two call shapes:
- `select(scope, selector) -> [atom]` — the set / equivalence class;
- `selectUnique(scope, selector) -> atom` IFF exactly one matches, else **FAIL
  LOUD** (log the N matches, return invalid). **No pick-first.**
This REMOVES the one-off's `ring.atomIndices.front()` anchor fallback (the latent
trap). The anchors become: PHE/TYR/HIS = `selectUnique(ring, {C, Gamma})`;
TRP-benzene = `selectUnique(benzene-ring, {C, Delta})` — unique *because the scope
is the ring* (CD1 is in the pyrrole ring, outside the benzene scope). Scoping is
what makes the typed selector unique without positional disambiguation. Day-1:
the scopes (`QtResidue.atomIndices`, `ring.atomIndices`) already exist; the
contract is the typed-filter-over-scope + the uniqueness semantics (a per-residue
`locant -> atoms` multimap is an optional convenience, not required).

**#4 — T1: bring it along, flagged unverified (Jessica).** Emit the FULL tensor
(T0 / T1 / T2) for every target; **FLAG T1** (and any cross-frame T2 component
where relevant) as `unverified` in the manifest/column metadata — the
Cartesian-NPY-vs-m-basis-H5 / 1e-vs-1o reconciliation is an open fileformat-schema
issue. Do NOT drop it (preserve-everything; `feedback_t2_sacred`; collect the
signal, flag don't discard).

**#5 — PBC: None (the protein is whole UPSTREAM).** The extractor cannot compute
anything without first putting the protein back together across the box; it does,
via the ported `pbc_whole` (`feedback_pbc_verbatim`, done in the extractor). So
the H5 positions arrive **whole**, every substrate query is protein-internal
(~≤10 Å, never solvent), and the substrate does **NO minimum-image**:
`PbcMode = None`. **Drop `PbcCellSeries`** from codex's proposal (8 → 7 files);
there is **no reader-side `pbc_whole.h`** (that's the extractor's job). Sole
obligation: a **load-time sanity assert** — fail loud if a frame's protein is
wrapped. The doc that kept misleading review agents (the "substrate MUST reuse
pbc_whole.h" mandate in `spec/substrate_conventions_2026-05-30.md`) is corrected
at the root.

**#6 — APBS units: pin from the producer.** `QtApbsGroup.h` says V/Å (E-field) /
V/Å² (EFG); `substrate_conventions` says e/Å² / e/Å³. Units cancel in a
correlation fit, but the manifest must label the true one — pin it from the
producer's APBS writer / H5 attrs before any manifest becomes an oracle artifact.

**#7 — record the frame per emitted relationship.** Charge multipoles are
atom-local; the new T2 work is lab-frame-equivariant. Each emitted relationship
records WHICH frame its tensor lives in (manifest/column metadata) — no consumer
guesses.

**#2 — input CLI: keep the existing flags. No JSON, no config framework.**
The input is the CLI we already have: `--run <trajectory>` (open) + `--out <dir>`
(write) + `--case <relationship>` (what to generate). The 9 relationships are new
`--case` values; one optional `--charge-source {ff14sb|aimnet2|mopac}` for the
multipole cases (or fold it into the case name). "Fail loud if a case needs
absent data" is the EXISTING behavior, not a new ValidateScenario system.
NOTE — do NOT conflate input with output: the per-relationship **output schema**
is the carrier above (which columns/NPYs each relationship emits — that is the
output concern); the **input** is trivial (open / which / where). The JSON-spec /
typed-params framework was over-reach; dropped.

## Layer 0 — data catalog + indexes (built DAY 1)

Everything the body holds sits in a typed **catalog**, one entry per array
(H5 dataset / NPY), so the **access language covers all of them uniformly** —
adding a datum is adding a catalog entry, not new grammar. (Mirrors
`python/nmr_extract/_catalog.py`'s one-`ArraySpec`-per-NPY discipline, on the
C++ body.) Components that "traverse every NPY" register here; the verbs in
Layer 1 then address any of them the same way.

```
ArrayId                                   # typed id per dataset/NPY
ArraySpec { id, rank, axes(atom?,frame?,slot?,comp?), dtype, accessor }
value(id, atom, frame[, slot][, comp]) -> scalar | vec3 | mat3   # one uniform addressing verb
```

Storage (resident, loaded once):
- `positions`    : `[atom][frame][3]`  **atom-major** (an atom's trajectory is
  contiguous — the inner-loop / Δ axis).
- `kernels[kind]`: `[atom][frame][tensor]` (ring / mc / efg / efield / charge … the
  H5 time-series; atom-major).
- `ring_nbhd`    : `[atom][frame][slot]{dist, rho, z, inplane, ring_id}`.
- `dft`          : sparse `(atom, orig_frame) -> {tensor3x3, orca_coord}`.
- `topology`     : static; atoms / bonds / rings / residues + reverse indices.

Indexes **built day 1** (at load, not lazily bolted on — scoping must be
O(log n) from frame 0 so the central loop never pays build cost mid-traversal,
and any frame, forward / back / random, is uniform-cost):
- `frame_map`: row↔orig, `dft_rows[]`.
- topology reverse indices: bonds-by-atom, rings-by-atom, ring-members (exist).
- **per-frame spatial KD-trees, ONE PER SOURCE-CLOUD KIND**:
  `atoms`, `bond-midpoints`, `ring-centers`, `charge-sites` →
  `trees[cloud_kind][frame]`. ~751 frames × few clouds × ~850 pts each is small
  (tens of KB/tree); build the set up front. Fallback if memory ever bites:
  lazy-but-cached; **default is day-1**.

Scoping is the point: a scenario picks a cloud + a cutoff and the tree answers a
**parametric radius query** — change the cutoff, no rebuild.

## Layer 1 — primitive access verbs (little APIs; pure; return views)

```
at(atom, frame)                       -> AtomState     # pos + kernel handles at (atom,frame)
window(atom, frame, w)                -> FrameSpan      # [frame-w .. frame]; ± by index arithmetic
neighbors(atom, frame, cloud, cutoff) -> [SourceRef]    # KD-tree radius query; cutoff REQUIRED, recorded
value(array_id, atom, frame, ...)     -> ...            # uniform catalog access (Layer 0)
```
- "Forward/back in time" = `window` / `±w` on the frame axis (atom-major ⇒ contiguous).
- "Move in space / scope" = `neighbors` with a parametric `cutoff` on the day-1 trees.
- No hidden cutoffs: `cutoff` is required at the call site and recorded (conventions).

## Layer 2 — relationship combinators (little-API lambdas)

```
Stratum        : (Body)                                      -> [atom]
LocalFrameFn   : (Body, atom, frame)                         -> LocalFrame
SourceSelector : (Body, atom, frame, cutoff)                 -> [SourceRef]   # rings | bonds | charges | ...
WindowFn       : (Body, atom, frame, w)                      -> FrameSpan     # time history (Δ / rate-of-change)
Attacher       : (Body, atom, frame, SourceRef, LocalFrame)  -> Values        # disp_local, kernel, n̂, identity
Reducer        : (SourceSet)                                 -> Record        # kernel form / sum-pool feature / Δ
TargetFn       : (Body, atom, frame)                         -> DftTarget
```
A relationship is a **named bundle**, defined in code:
```
Relationship { name, stratum, frame_fn, selectors[], attachers[], window_fn?, reducer, target_fn }
```
- **Heterogeneous sources** = `selectors = [rings, bonds, charges]`. No special
  case — just more lambdas. (The breadth-via-mechanisms need, structurally.)
- **Differencing / rate-of-change** = supply `window_fn` + a reducer that takes Δ
  over the window. A combinator, not a rewrite. (The identifiability test rides in.)

## Layer 3 — the engine (one pure, order-free loop) + sink

```
for (atom, frame) in index_space:            # order is FREE (pure transforms) → choose for cache
    lf   = frame_fn(body, atom, frame)
    src  = flatten( sel(body, atom, frame, cut) for sel in selectors )
    vals = [ attach(body, atom, frame, s, lf) for attach in attachers for s in src ]
    rec  = reducer( SourceSet(src, vals, window_fn?(body, atom, frame, w)) )
    sink.fold(rec)
```
Schedule is a **cache** choice, not a correctness one (referential transparency):
per-frame-outer reuses a frame's KD-trees across atoms; per-atom-outer keeps an
atom's Δ-window contiguous. Both O(1)-amortized given day-1 trees + atom-major
storage — pick per profile.

Sink = streaming **fold** (CSV now; in-memory / pybind-exposed buffer later). No
materialized cross-product.

## Performance notes
- atom-major storage: inner-loop time-window / Δ is a contiguous span; per-frame
  trees are built once (day 1), queried many times.
- type erasure: the innermost `Attacher` runs N_src × N_atom × N_frame; if
  `std::function` ever bites, **template the attachers** (selectors / frame /
  reducer are nowhere near hot enough to care). Batch process → fine to start
  type-erased.

## The functional language, grounded in the one-off

Read off what the one-off actually does (cited), so each verb generalizes a real
operation rather than inventing vocabulary. Two sides share the SAME core nouns
(group / window / value / classify) so the language folds across the C++/Python
line.

### Producer side (C++) — primitive verbs (pure, return views)
| verb | signature | generalizes (one-off) |
|---|---|---|
| `atomsWhere` | (pred) -> [atom] | stratum loops `IsAromaticRingHydrogen` / `IsBackboneAmideHydrogen` (RingCurrent/McConnell `extract`) |
| `pos` | (atom, frame) -> vec3 | `conf.atomPosition` (atom-major) |
| `window` | (atom, frame, w) -> FrameSpan | **NEW** — the time axis (one-off had no Δ); contiguous via atom-major |
| `value` | (arrayId, atom, frame[,slot][,comp]) | `rn->at`, `bs->at`, `mc->at` → one uniform catalog read |
| `near` | (cloud, atom, frame, cutoff) -> [SourceRef] | `FrameSpatialIndex.Within(pos, cutoff)` (bonds) → per-cloud |
| `slots` | (slotList, atom, frame) -> [SourceRef] | the `ring_nbhd` slot walk `rn->ringIndexAt`/`rn->at` (frozen membership) |
| `atomOf` | (residue, locant[,branch]) -> atom | `typedRingAnchor`'s Locant scan + `QtResidue` N/CA/C cache → typed O(1) |
| `ringsOf` | (atom) -> [ring] | `ownAromaticRing` / `ringMembershipsForAtom` |
| `ringGeom` | (ring, frame) -> {center,normal,radius} | `RingGeometryAt` + the canonical-normal flip (now a body op) |

### Producer side — relationship combinators (the little-API lambdas)
| combinator | generalizes (one-off) |
|---|---|
| `Stratum = atomsWhere(pred)` | the typed stratum filters |
| `LocalFrameFn` | ring: `ringsOf→ringGeom→canon-normal→atomOf(res, Γ\|Δ)→BuildAromaticHFrame`; HN: N/H/CA/Cprev `→BuildHNFrame` |
| `SourceSelector` (2 backends) | `slots(ring_nbhd)` (ring) **or** `near(bond-midpoints, cutoff)` (mc). Unify both. |
| `SourceClassifier` | `is_self_or_bonded` (ownRings/ownAtoms overlap) + the type key (ring_type / bond_category) |
| `Attacher` (geometry) | `disp_local=toLocal(center−pos)`, `r`, `cosθ`, `(3cos²−1)/r³`, `source_normal_local`, `bond_axis_local` |
| `Attacher` (identity) | `QtRing`/`QtBond` virtuals (intensity, NitrogenCount, JBoffset, category, order, elems) |
| `Attacher` (kernel) | `value(bare_kernel_TS, atom, frame)` — the producer cross-check |
| `WindowFn` | **NEW** — `window(atom,frame,w)` for Δ-features (the rate-of-change test) |
| `Reducer` | `sumAll` / `sumValid`(filter !self) / per-key sums / counts; or pass-through (un-summed rows); or Δ over window |
| `TargetFn` | `BuildTarget`: `dft(atom,orig)→raw3×3→DecomposeLibrary(T0/T1/T2)→toLocal(total_local)` |

### Consumer side (Python) — the SAME nouns, named once (replaces the pile)
| verb | generalizes (the script pile) |
|---|---|
| `by(keys)` | `groupby(atom,frame)` / `(atom,ring)` (look02, equiv_t2) |
| `demean_within(key, col)` | within-atom de-mean / fixed-effects baseline strip (look01, credibility, equiv_t2) |
| `windowed_diff(series, w)` | finite difference (diag_differencing) — Δ |
| `autocorr1(series)` | lag-1 smoothness gate (diag_differencing) |
| `fit(x,y) / leave_out(group)` | corr / slope / R² / leave-atoms-out (credibility, credibility2) |
| `libT2(M3×3) -> vec5` | the library-basis spherical projection (equiv_t2; matches emitted to 4.9e-8) |
| `sumpool(sources, group, fn)` | scatter-add pooling + per-group baseline (sumpool_*, equiv_t2) |

The pile becomes `load(substrate) |> by |> demean_within |> {fit | sumpool}`. The
shared nouns (`by`, `window`/`windowed_diff`, `value`, classify) mean the access
language is the same grammar on both sides of the CLI line — which is the point.

### Iterated closures (keep this — it's the heart, not a struct-of-functions)

The verbs are not free functions you pass args to each call; they are **closures
that capture the body + their config (currying)**, and the engine **iterates**
them. Partial application builds the specialized verb; iteration applies it; the
fold is a state-carrying closure. That is what makes "a set of lambdas each with
a little API" actually *compose* — lose it and you are back to the one-off's
hand-inlined procedural loop.

Build (each line is a closure with config + body captured once):
```
sel    = near(BondMidpoints, cutoff=8.0)     # closure (atom,frame) -> [SourceRef]   (cloud+cutoff baked in)
       | slots(ring_nbhd)                     # ...or the precomputed-slot backend
frameF = aromaticHFrame                       # closure (atom,frame) -> LocalFrame
attach = [ disp_in(frameF), normalIn(frameF), kernelOf(BS) ]   # closures (atom,frame,src) -> Values
classify = selfBonded(ringsOf) ⊕ ringType     # closure (atom,frame,src) -> {valid?, key}
reduce = sumByKey(classify)                    # STATE-CARRYING closure (fold)
```
Iterate (lazy, two levels — closures producing closures; pure ⇒ order-free):
```
fold( sink,
      map( λ(atom,frame): reduce( map( attach, sel(atom,frame) ) ),   # inner: over sources
           index_space ) )                                            # outer: over (atom,frame)
```
- **Currying** is the composition mechanism: `near(cloud,cut)` returns the
  configured `(atom,frame)->[src]` closure; `disp_in(frameF)` captures the frame
  closure; a relationship is just these curried closures bound together. No
  args threaded through the loop, no body passed everywhere.
- **Iteration** is map over the index space + an inner map over sources + a fold.
  The one-off *is* this with every closure hand-inlined and un-curried; the
  language curries the verbs and lets the engine iterate — same computation,
  now composable and reorderable.
- **C++ reality:** closures capture `const Body&` + config; cool levels
  (relationship / frame / selector) may be `std::function`; only the innermost
  per-source closure is templated *if* profiling bites. Python: plain closures /
  `functools.partial`, identical shape. Same idea both sides of the line.
- **WindowFn / Δ** is just another captured closure (`window(·,·,w)`) folded in —
  iterated closures absorb the rate-of-change test without new structure.

### Additional indexing tools to posit (what the verbs need, built day 1)
1. **Per-cloud per-frame KD-trees**: `atoms` / `bond-midpoints` / `ring-centers`
   / `charge-sites`. Generalizes the single bond-midpoint `FrameSpatialIndex`.
   Powers `near(cloud, …)` with parametric cutoff. Built day 1.
2. **Typed-atom index** `(residue, locant[, branch, element]) -> atom`.
   Generalizes `QtResidue`'s N/CA/C/O/H/HA/CB cache to ALL locants. Powers
   `atomOf` (frame anchors, typed source/target selection) with no name-scan and
   no positional guess — the identity-from-chemistry index.
3. **Own-ring / own-atom-set per aromatic-H** (for `is_self_or_bonded`): the
   one-off rebuilt this per run; make it a day-1 index off `ringMembershipsForAtom`.
4. **Slot-list backend** (`ring_nbhd`): already an index in the H5; just name it
   as the `slots()` source backend (no build).
5. **Atom-major layout IS the time-window index** — `window(atom,frame,w)` is a
   contiguous span; no new structure, note it.
6. (consumer) group-offset index `(atom,frame)`/`(atom,ring)`: pandas now; a
   body-side view later only if the fitter binds in-process.

## Python interface (RESOLVED — it's a batch file, don't hook it up)

Mental model: a **batch file**. A shell script of CLI invocations, each writing
its log to a file and its results to disk. No UDP, no sockets, no live channel —
that's GUI thinking. The interface is a **process + file contract**.
- **IN**: a JSON **scenario spec** — selects a *named* relationship + params
  (calcset, strata, cutoffs, window, out). Python (or the script) composes it;
  the CLI reads it. (Declarative format ⇒ discuss before finalizing,
  `feedback_data_format_discuss_first`.)
- **RUN**: the current **synchronous** `main` (load → traverse → commit → rc).
  No event loop / worker needed. Fired from the batch file:
  `h5reader_extract --scenario spec.json --out dir > dir/run.log 2>&1`.
- **LOG-BACK**: the per-run **`.log` file** (the `StructuredLogger` already
  writes stderr; the batch redirect captures it). UDP 9997 is the *GUI's*
  live-tail channel — **irrelevant to batch; do not design around it.**
- **OUT**: atomic `QSaveFile` substrate + a `manifest.json` ("the end bit").
  `QSaveFile` commits **only on success** ⇒ output-exists-iff-success.
- **DONE/OK**: process **exit code** + manifest present. (`subprocess.wait()` /
  the script's `$?` IS the process protocol — NOT mid-run file-polling.)

So the four original asks collapse to almost nothing: **log-back** = the `.log`
file; **clear-exit** = exit code + atomic output + manifest; **commands-in** and
**events-as-IPC** = **dropped** (a batch job carries everything in its spec up
front and runs to completion). Python = compose spec → run (or let the batch
file run) → read manifest + substrate; **never opens `trajectory.h5`**.

Codex (read-only research, task `byzrc0uff`) recommended a heavier control plane
— a `QLocalSocket` command/event channel, a v1 event schema, cooperative
cancellation, and reshaping `main` into `QtConcurrent`+`QPromise`+`QFutureWatcher`.
That is justified ONLY by **mid-run** needs (cancel a long run; stream
authoritative live events). A batch file needs none of it, so **not adopted
now**. Crude-correct cancel = kill the process (atomic commit ⇒ no corruption)
until a run is long enough to want cooperative cancel.

**Kept from codex regardless:** a versioned JSON schema for the manifest + final
log line; the rc classes (2 bad-args, 1 load, 3 sink, 4 commit, 130 cancel); the
Qt 6.4 notes (`QPromise`/`QLocalSocket` present; avoid `QNativeIpcKey`,
`Q_APPLICATION_STATIC`); and the **documented upgrade path** — codex's
QLocalSocket+QtConcurrent design (same schema/run-id/cancel states) migrates in
verbatim if/when a run needs live control.

**End-run anti-patterns (adopt wholesale, from codex):** Python re-reading
H5/NPY; scraping stdout/stderr for state; polling output files mid-run for
completion; UDP-only authoritative events; a second logger; duplicating model
state in Python; per-row Qt signals; partial writes to final paths (use
`QSaveFile`); `processEvents()` in hot loops; `QThread::terminate()` as normal
cancel.

## Validation
- each lambda unit-tested against its little API in isolation;
- the composed **ring-current** relationship reproduces the one-off **ORACLE**:
  leave-atoms-out k≈21 ppm·Å³, R²=0.62 (coupled); equivariant T2 R²=0.44; |T2|
  r=0.75; basis check 4.9e-8; DFT frame rotation ~1e-4°. Match ⇒ faithful rebuild.

## Constraints carried (don't move)
2 proteins have trajectory DFT — precious; **T2 only from DFT** (NMR is
isotropic); **fleet = NMR/isotropic**, separate. Breadth via **mechanisms within
a protein, not many proteins**. Identifiability limited — differencing stays
**open** (technique ≠ identity). Instantaneous model; correlate-not-match.
**Reader owns H5** (Python never opens it). Build/verify in the lead session;
additive to the reader; GUI untouched; experimental one-shot branch.

## Agent briefs (the serious version) + the gate rule (Jessica, 2026-05-31)

A stub / implementation agent for this surface must be GIVEN (not left to
discover), or it cannot do the work seriously:
- what is in the **OBJECT MODEL** (`OBJECT_MODEL.md` + the reader's typed model:
  `QtProtein`/`QtConformation`/`QtFrame`, the `QtAtom` hierarchy + predicates,
  the `QtRing` hierarchy + virtuals, typed residues/bonds, `Mat3`+`SphericalTensor`);
- what is in the **NPYs / H5 datasets, described in terms of the reader OBJECTS**
  they become (not raw array names);
- the **MATHS of what we are doing and WHY, both levels** — scalar kernels up to
  equivariant T2 — with the physics rationale;
- a directive to **FOLLOW the iterated-closure / lambda model** over the
  fully-loaded, in-memory-indexed trajectory — NOT to tear it apart or
  second-guess it. State why it is the right shape for our maths on a resident
  indexed trajectory, and have the agent follow it out.

Process (hard rule): **propose the brief and GATE on confirmation before
firing.** Never launch an agent while the brief is still being written — it
burns ~100k tokens for nothing. (Learned the expensive way, 2026-05-31.)
```
