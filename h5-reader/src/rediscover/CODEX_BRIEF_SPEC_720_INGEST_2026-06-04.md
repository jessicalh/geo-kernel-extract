# Codex brief — SPEC the 720-WT static-pose ingest (design doc only, no build)

> **Historical session brief — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Status: **DRAFT — pending lead plan-vet, not yet fired.**

You own the grind; the lead vets + judges + owns ALL git. This is a **spec-authoring** pass, not a
build. We are stashing our understanding of the 720-WT statics work as a precise design document the
lead will vet and then shelve ~2 days. Your deliverable is ONE markdown spec. Write no code, run
nothing, touch no substrate. Honest collaborator stakes: a clean spec now is what lets us execute the
720 ingest in one disciplined loop later instead of re-deriving it cold.

## CONTEXT (read first, in order)
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/NEXT_SESSION_PROMPT.md` (STEP 3 B-path),
  `NOW.md`, top of `STATE.md` — what the 720 pilot is FOR (the only between-axis instrument; 1P9J is
  a within instrument; the seed of the unified stats engine).
- `DESIGN.md` + `GUIDANCE.md` (this dir) — the rediscover class model (`RunData`/`RunLoader`/`FrameMap`/
  `DftFrameSet`/`LocalFrameBasis`/`SpatialIndexSet`/`FrameSpatialIndex`/`TypedAtomIndex`/
  `RingGeometryCache`/`ResidentIndexes`), and the "model is the spine, emit C++-side, Python only
  consumes" boundary.
- `PerAtomSubstrate.h` + skim `PerAtomSubstrate.cpp`, `Catalog.{h,cpp}`, `RunData.{h,cpp}`,
  `ExtractionSupport.{h,cpp}`, `main_extract.cpp` — the current per-(atom,frame) emit + the trajectory
  load path you are mirroring for statics (minus the time axis).
- The reader's NPY/topology loaders (`src/io/QtNpyReader`, `src/io/QtTopologySidecar`, `src/io/QtProteinLoader`,
  `src/io/DftShieldingLoader`/`OrcaShieldingParser`) — how the model is populated from the producer's
  NPY surface + raw DFT tensors.
- `PARTITION_FILTER_DESIGN.md`, `NODE_STORE_CONTRACT_2026-06-02.md` — the substrate column contract.
- Memories: `project_unified_stats_engine`, `feedback_all_statistics_minimize_python_15gb`,
  `feedback_model_is_spine`, `feedback_no_file_discovery`, `feedback_static_workaround_not_producer_redo`.

## THE INPUT — a GIVEN; design against it, do NOT hunt for it
The canonical input is **the curated 720 mutant set** (the 725 in `/shared/2026Thesis/consolidated/` is
OUT OF DATE — has duplicates; do not target it). Its exact path is **TBD** (the set and/or its NPYs are
likely in `learn/`, currently recovering from the spinner). **Do NOT block on the path, and do NOT assume
co-located NPYs.**

**Assume as a GIVEN for this plan:** the extractor will be **re-run over all 720 (curated), MOPAC ON**, so
the full per-protein producer output — the canonical NPY/H5 surface (kernel shadows; APBS + AIMNet2 + MOPAC
fields; the raw DFT dia/para/total tensors = the **WT absolute σ** target; positions; topology) — **exists
per protein, keyed by protein ID.** The ingest **CONSUMES that per-protein producer output**, exactly as
the trajectory path consumes the trajectory H5, minus the time axis. The extractor is run **UNMODIFIED**
(SACRED). Reference inputs by the documented per-protein-ID convention; log-and-stop on any missing expected
file — no glob, no try-and-fail (`feedback_no_file_discovery`). The exact filesystem path is a later
execution fill (when `learn/` recovers); the spec is path-agnostic.

## THE FORECLOSURE — the #1 thing this spec must rule out (the lead's hard line)
The spec must **explicitly and structurally foreclose** two failure modes — name them and design them out,
do not merely discourage:
1. **No terabyte of Python.** No bulk / per-source / pairwise materialization to Python; no boa-constrictor.
   The emit is the lean keyed digest (< 15 G total) and ALL reductions happen C++-side. (The prior 57 GB
   pairwise/occupancy dump is the named anti-pattern; a 720× version is a terabyte — forbidden.)
2. **No second protein model in Python — not even a secondary aggregate one.** The spatial work (KD-trees,
   neighbour finding, geometry, ring/bond discovery) and the protein/topology model live in the typed C++
   model + its indexes, FULL STOP. Python never holds proteins/atoms/residues as a model, never does spatial
   work, never builds an "aggregate" model to stitch the 720 together. Python's ONLY role is the fitter over
   a lean tabular substrate.

**Pre-empt the rationalization.** If it ever feels like "Python just needs the model/spatial structure to do
X — I had to do it this way," that is the SIGNAL to extend the C++ emit or do the reduction C++-side, NEVER
to rebuild the model or the spatial work in Python. The spec states this as law and shows the design honors
it on BOTH memory-strategy branches (the option-(b) "buffers of the things we care about" are C++-side
accumulators, not a Python aggregate).

## WHAT THE SPEC MUST DESIGN
A `StaticPoseConformation` ingest path in the rediscover extractor that, over the curated 720:
1. loads each protein's per-protein producer output into the typed C++ model (reusing `QtProtein`/topology +
   the resident indexes — `SpatialIndexSet`/`FrameSpatialIndex`/`TypedAtomIndex`/`RingGeometryCache`), builds
   the local frames — a static pose is effectively one frame, so most per-(atom,frame) machinery carries over;
2. accumulates the rediscover substrate the fitter consumes — the source neighbourhoods (displacement vectors,
   local frames) computed by the model's spatial tools from positions; the kernel shadows read from the
   producer output; the DFT target from the raw dia/para/total tensors;
3. emits the **per-(protein, atom)** substrate — the same column contract the fitter consumes — under the
   memory strategy the maths dictates (fork below);
4. uses the typed C++ model + indexes for ALL accumulation; **no second Python model**; **lean disk < 15 G
   total output**; generous RAM (128 GB, swap OK). Python only consumes the lean emit.

## THE DESIGN FORKS — lay out, recommend where you can, FLAG for the lead (do NOT bake)
**(1) Kernel sourcing — RESOLVED (given).** The per-protein producer output exists (the 720 re-run with
MOPAC). The ingest consumes it; no `nmr_extract` re-implementation, no kernel recompute reader-side. Adapt
the per-protein NPY/H5 surface into the typed model as the trajectory path does. Be concrete about what the
static path sources from where: positions, topology, the kernel shadows, the DFT target (raw dia/para/total).
**(2) Memory / accumulation strategy — OPEN, depends on the maths (do NOT bake; the lead is explicit it is
undecided — it follows the between-axis statistics method, not yet settled).** Lay out BOTH, design so
neither is foreclosed, and recommend a structure that supports both with the binding decision deferred:
- **(a) all-720-resident → emit once:** hold the whole curated set in the C++ model + cross-protein indexes,
  then emit the lean digest at the end. Enables between-protein statistics that need all proteins co-resident.
  Costs RAM (acceptable here — 128 GB, swap OK).
- **(b) sequential, accumulate running buffers:** stream protein-by-protein, accumulating ONLY the reduced
  quantities the between-axis maths needs ("buffers of the things we care about" — the boa-constrictor
  avoidance), emit at the end. Lighter RAM.

## THE ANTI-REGRESSION GATE
Define an **oracle-parity gate**: a known 1P9J WT static pose run through the static path must reproduce the
corresponding 1P9J trajectory-frame substrate values (or a precisely-stated equivalence, since a static pose
≠ a sampled frame). State exactly what must match and to what tolerance, and what may differ and why. The
gate is how we trust the static path before trusting 720 numbers.

## ALSO COVER
- The emit column set for statics: which per-(atom,frame) columns carry over, which are trajectory-only
  (Welford/frame-modulation) and drop, which are new (per-protein identity / provenance). The lean menu the
  fitter needs — surface the format for lead vet (`feedback_data_format_discuss_first`), don't invent schema.
- A per-protein run manifest (status / atom counts); support-flag thin element strata.
- Disk budget arithmetic (rows × cols × bytes) showing < 15 G, under both memory strategies.
- BUILT vs DESIGNED vs OPEN vs CAVEATED throughout; an "Open / undecided" running list (the memory-strategy
  fork lives there, tied to the maths).

## After you: an opus agent adversarially reviews this spec
An opus agent will sharpen the memory-strategy fork (keep it honestly open + tied to the maths, not
prematurely resolved), check the consume-the-given-producer-output path is concrete, pressure-test the
oracle gate, and guard the model-is-spine / lean-disk / no-second-Python-model constraints. Write it to survive that.

## HARD RULES
- **DOCS ONLY.** Write exactly ONE file:
  `/shared/2026Thesis/nmr-shielding/h5-reader/src/rediscover/SPEC_720_STATIC_INGEST_2026-06-04.md`.
  Change NO code, run NO build / extraction / 720 / ORCA / fit, touch NO substrate or H5, run NO `git`.
  Everything else is READ-ONLY.
- **No file discovery in the design** — documented per-protein-ID convention only; log-and-stop on missing.
- **Encode the constraints as law in the spec — THE FORECLOSURE above is the headline:** model-is-spine /
  **no second protein model in Python (not even a secondary aggregate)** / **no spatial work in Python** /
  **no terabyte Python dump** / Python only consumes the lean emitted substrate, reductions C++-side; lean
  DISK < 15 G, generous RAM; extractor SACRED (run unmodified); frozen `get_C`; never open `trajectory.h5` in Python.
- **The memory-strategy fork stays OPEN** (maths-dependent) — recommend a non-foreclosing structure, do not pick.
- Branch `h5-reader-pysr-spike` — never merge/switch/rebase/PR. **Lead owns ALL git.**
- Truthful, no overclaiming; cite `path:line` for load-bearing claims about existing code.
