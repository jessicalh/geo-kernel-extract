# Post-hoc DFT → `.LGS` registration bridge — scope & design

**Date:** 2026-06-05
**Status:** SCOPE / DESIGN — and **SINCE EXECUTED (2026-06-05):** acting on this design, the dense
LGS was re-registered 660→751 DFT frames and the lean LGS stashed (`.LGS.back`); the bridge is now
operational (see the spine doc CURRENT STATE). The prose below is the original design pass.
**Scope licence:** reader stays read-only; T2 preserved end-to-end; the
model is the spine; one registration serves BOTH consumers.

This note answers the lead's question: register the now-complete 1P9J
ORCA campaign in the `.LGS` so the DFT shows up in BOTH the rediscovery
extractor AND the h5-reader display (`orca_dft:*`), from ONE shared
registration. The headline up front, then the evidence, then the flow.

---

## Headline

**It is a read-only reader-wiring job, NOT a producer-extraction step.**
The reader already parses a completed ORCA campaign straight into
`DftShieldingStore` at load — `meta.json → files.out_primary → parse the
.out → validate → resident frame` — with no H5 write and no re-extraction
(`DftShieldingLoader::LoadAndValidate`, `io/DftShieldingLoader.cpp:45`).
`orca_dft` is NOT a dataset inside `trajectory.h5`; it is a frame-local
source the reader reads from disk on demand. Nothing needs to be written
into the H5, and the extractor is not touched.

The single registration that already serves both consumers is the
`.LGS`'s `dft.frames[]` array (the typed `frame_index → meta_json` map).
The only thing stale is the *contents* of that array: the campaign grew
from 660 to 751 completed frames after the current `.LGS` was written.
"Making up the difference" is re-running the existing registrar
(`tools/lgs_write.py`) over the completed jobs — it is already idempotent
and already filters to completed jobs. There is a second, operational
gap (two `.LGS` in one calcset dir) covered in §6.

---

## 1. Why `orca_dft:*` reads empty — the exact gap

### 1a. The read path and what it expects

The reader's DFT chain, end to end:

1. `CalcsetManifest::Load` parses the `.LGS`; the optional `dft` block
   yields `dft.frames[]`, each a `DftFrame { frame_index, meta_json_abspath }`
   (`io/CalcsetManifest.cpp:459-526`; type at `io/CalcsetManifest.h:48`).
2. `ReaderMainWindow` builds the store **only if** `manifest.dft` is
   present, straight from `dft.frames` — no dir scan, no name parsing
   (`app/ReaderMainWindow.cpp:467-469`).
3. `DftShieldingStore` copies `frame_index → meta_json_abspath` into a
   hash; `.out` files parse lazily on `requestFrame()`
   (`model/DftShieldingStore.cpp:18-34`, `.h:67`).
4. Per frame, `DftShieldingLoader::LoadAndValidate` opens the meta.json,
   reads `files.out_primary` (basename, resolved against the job dir —
   NOT a glob), parses the `.out`, and **rejects** the frame on any of:
   cannot-open, no `out_primary`, parser hole, `atoms.size() !=
   protein->atomCount()`, or `total != dia+para` to 0.1 ppm
   (`io/DftShieldingLoader.cpp:45-111`).
5. The catalogue advertises `orca_dft:{total,diamagnetic,paramagnetic}`
   unconditionally — the descriptor exists whenever `storagePath`
   starts with `orca_` (`model/TrajectorySignalCatalog.cpp:250-251`,
   registered at `:1037-1061`). So the *channels* are always present in
   the picker; their *samples* are empty whenever the store has no
   resident frame for the current original index.

**So `orca_dft:*` reads empty in exactly two situations:** (a) the store
was never built because `manifest.dft` was absent/`nullopt`, or the whole
load failed; or (b) the store was built but every `requestFrame` resolved
absent (a validation reject, or no job for that original index).

### 1b. What the "dense LGS" actually contains (verified)

The dense LGS named in the brief —
`…/1p9j-calibration-with-dft/1p9j-calibration-dense-mopac-live-orca.LGS`
— is **well-formed and its `dft.frames[]` is populated**, not empty:

- `schema_version: 1`, `kind: trajectory`, `dataset_id:
  1p9j-calibration-dense-mopac-live-orca`, `protein_id: 1P9J_5801`.
- `trajectory` points at the dense all-frame MOPAC extraction
  (`extract-full-mopac-allframes-20260531/…`), `frame_index_basis:
  trr_frame_index`, `frame_dt_ps: 10.0`.
- `dft.frame_stride: { first: 0, last: 1318, step: 2 }`,
  `campaign_target_frames: 751`.
- `dft.frames[]` has **660 entries**, even `frame_index` 0 → 1318, each
  `meta_json` an absolute path into `/shared/2026Thesis/1p9j-orcas/jobs/…`.
- The file closes cleanly with a valid `metadata` block
  (`generated_at_utc: 2026-06-02T11:43:12Z`, `lgs_writer: lgs-tools 0.1.1`).

The pointed-at meta.json files have the exact shape the loader needs:
`files.out_primary` (a basename in the job dir), `frame_index`,
`frame_ps`, `orca_exit_code: 0`, `atom_count: 846` (e.g.
`…_f000000_t0.0_meta.json`). The matching `.out` exists and has **846
`Nucleus` records inside the `CHEMICAL SHIELDINGS` section** — equal to
the dense extraction's topology atom count (846), so the atom-count gate
passes and the `total == dia + para` identity holds (ORCA prints all
three tensors). A frame loaded from this LGS would validate and become
resident.

**Therefore the empty-`orca_dft` symptom is NOT "the dense LGS fails to
fire the loader on its own frames."** Loaded by its explicit `.LGS` path,
the dense LGS populates the store with 660 valid frames. The real gaps
are two, and neither is a parser/loader defect:

- **Gap A — staleness (the post-hoc problem in the brief title).** The
  campaign is now **complete: 751 jobs, even frames 0 → 1500**
  (`/shared/2026Thesis/1p9j-orcas/jobs/` has 751 dirs; the last is
  `…_f001500_t15000.0`; `status/jobs.csv` is header + 751 rows). The
  dense LGS was written at frame 1318 and so omits 91 completed frames
  (1320 → 1500). Every consumer that reads this LGS silently under-covers
  the back third of the trajectory. This is the structural post-hoc
  issue: MD → extract → write `.LGS` → ORCA keeps finishing → the `.LGS`
  is now behind the campaign.

- **Gap B — directory ambiguity (the likely "reads empty" trigger).**
  The calcset dir now holds **two** `.LGS` files
  (`1p9j-calibration-with-dft.LGS` — the older, sparse-`extract/`,
  500-frame one — and the dense one above). `CalcsetManifest::Load`,
  when handed a **directory**, requires exactly one `*.LGS` and treats
  zero-or-many as a hard error (`io/CalcsetManifest.cpp:166-177`,
  spec `CALCSET_MANIFEST.md:46-48,184-187`). So "open the calcset
  directory" now fails the whole load → no protein, no store, and the
  `orca_dft:*` channels are present-but-empty. If instead the reader was
  pointed at the *old* `1p9j-calibration-with-dft.LGS` by file path, the
  store builds but DFT stops at frame 1002 (500 frames) — under-coverage,
  not emptiness. Either way the fix is registration hygiene, not loader
  code.

### 1c. The form the loader expects (restated, so the registration matches it)

`DftShieldingLoader` expects, per frame: a `meta_json` path that exists,
whose `files.out_primary` names an `.out` sitting next to it, whose
`CHEMICAL SHIELDINGS` section has one `Nucleus` block per topology atom,
with `Diamagnetic`/`Paramagnetic`/`Total shielding tensor` 3×3 matrices
each (`io/OrcaShieldingParser.cpp:46-126`). The campaign's meta.json +
`.out` already satisfy this. The registration's only job is to list the
right `frame_index → meta_json` pairs; it must not duplicate or transform
any per-frame content (spec `CALCSET_MANIFEST.md:165-178`).

---

## 2. Read-only? (the headline, evidenced)

**Yes — the reader parses the completed campaign into `DftShieldingStore`
on load; no H5 write, no re-extraction is required.**

- The store is a per-frame *source provider*: `requestFrame(orig)` parses
  one `.out`, makes it resident, emits `frameReady`; the next request
  releases the prior frame. Persistent history lives in the strip
  `ChannelBuffers`, not in the store, and never in the H5
  (`model/DftShieldingStore.h:18-38,88-126`; `.cpp:46-93`).
- The only file I/O is `QFile`-read of meta.json + `.out`
  (`io/DftShieldingLoader.cpp:49-75`). No `trajectory.h5` handle is
  opened for DFT; `orca_dft` is `SignalSourceKind::OrcaDftFrame`,
  explicitly NOT an H5 dense path (`TrajectorySignalCatalog.cpp:250`).
- This honours the standing rules: viewer never writes H5 / never
  re-extracts (`h5-reader/CLAUDE.md`, memory `feedback_qt_read_only_h5`),
  and Python never opens `trajectory.h5` (memory `no_parallel_h5`) — the
  reader is the one that reads it, and here it doesn't even need to.

The T2 surface is preserved throughout: the parser keeps the raw 3×3 AND
the `SphericalTensor` (T0/T1/T2) for total/dia/para; the store samples
either T0 (isotropic) or |T2| (anisotropy) per channel, never collapsing
the tensor (`model/DftShielding.h:20-53`,
`model/DftShieldingStore.cpp:73-86`). No producer extraction step exists
or is needed for the reader to show DFT.

(The contrast worth stating plainly: an extraction step would only be
needed if `orca_dft` were a *baked H5 dataset*. It is not — by design it
is a sidecar the reader resolves live. That design is what makes this
read-only.)

---

## 3. Shared registration — one shape, both consumers

There is **already exactly one registration**, and both consumers already
read it. No second scheme is needed or wanted.

**The single registration is `manifest.dft.frames[]`** — the typed
`frame_index → meta_json` map in the `.LGS` (spec
`CALCSET_MANIFEST.md:85-95,142-178`).

- **Rediscovery** pulls DFT by walking `manifest.dft->frames` and calling
  the SAME `DftShieldingLoader::LoadAndValidate(f.meta_json_abspath,
  protein)`, keying the result by `f.frame_index` into `DftFrameSet`
  (`rediscover/RunData.cpp:120-130`). `FrameMap::Build` then joins DFT to
  H5 rows purely by original (trr) index, after a hard `frame_index_basis
  == "trr_frame_index"` check — so a sparse/irregular even-DFT set over a
  dense stride-1 extraction maps correctly with no extra schema
  (`rediscover/RunData.cpp:46-74`; grounding doc
  `rediscover/MULTISTRIDE_LGS_GROUNDING.md` §1, §4). The substrate then
  carries `dft_present` + the raw/decomposed tensors per atom
  (`rediscover/ExtractionSupport.cpp:45-68,164-185`).
- **Reader** builds `DftShieldingStore` from the SAME `dft.frames`
  (`app/ReaderMainWindow.cpp:467-469`) and samples it for the strip
  chart, keyed by the SAME original index
  (`model/DftShieldingStore.cpp:36-86`).

Both go through `DftShieldingLoader::LoadAndValidate` and both key on
`frame_index == /trajectory/frames/original_index`. The registration is
method-agnostic about *who* reads it: it is just the index. So the design
rule is simply **"keep `dft.frames[]` complete and correct, in one
canonical `.LGS` per calcset"** — and both consumers light up from it.

What the registration must carry for both (already in v1, no schema
change):

| Field | Why both need it |
|-------|------------------|
| `dft.frames[].frame_index` | the join key to H5 `original_index` (reader resident lookup; rediscovery `FrameMap`) |
| `dft.frames[].meta_json` | the per-frame pointer both pass to `DftShieldingLoader` |
| `trajectory.frame_index_basis` | rediscovery hard-checks `trr_frame_index`; reader maps rows on it |
| `dft.method`, `dft.campaign_target_frames`, `dft.frame_stride` | display + glance-level cadence; coverage still comes from `frames[]` |

No per-consumer field, no second map, no DFT-to-extraction-row table is
required: the original-index join already spans dense extraction × sparse
DFT (grounding doc §5 confirms v1 is sufficient for this exact
dense/sparse case).

---

## 4. The flow: REGISTER → CHECK → MAKE UP THE DIFFERENCE

The post-hoc campaign keeps finishing frames after extraction. The flow
must pick up newly-complete frames on each load/run, idempotently, and
honour multi-stride (dense stride-1 MOPAC extraction + sparse stride-2
DFT in one `.LGS`).

### REGISTER (campaign → `.LGS`)

`tools/lgs_write.py` IS the registrar and already does the right thing:

- Walks `dft/jobs/*/` or an explicit `--dft-jobs-dir`, reads each
  `*_meta.json` for the typed `frame_index` (NOT the dir name), **filters
  to `orca_exit_code == 0`**, sorts by `frame_index`, and emits
  `dft.frames[]` plus a derived `frame_stride`
  (`tools/lgs_write.py:145-219,283-323`). The current dense LGS was made
  this way with `--dft-jobs-dir /shared/2026Thesis/1p9j-orcas/jobs`.
- Output is `<root>/<dataset_id>.LGS`; refuses to clobber without
  `--force` (`tools/lgs_write.py:444-453`).

This is the registration both consumers read. It is the multi-stride
registration already: `trajectory` points at the dense stride-1 MOPAC
extraction; `dft.frames[]` lists the sparse stride-2 even frames; the two
strides coexist in one `.LGS` and join by original index (no
`extraction_stride` field needed — grounding doc §5).

### CHECK (parsed/extracted vs available)

Two honest comparisons, both already computable, neither requires new
schema:

- **Coverage check (registration vs campaign):** the campaign target is
  `dft.campaign_target_frames` (751) and the planned cadence is
  `dft.frame_stride`; the *registered* set is `len(dft.frames[])` (660);
  the *available* set is the count of `orca_exit_code == 0` jobs on disk
  (751). `registered < available` ⇒ the `.LGS` is stale ⇒ re-register.
  `lgs_write.py` recomputes the available set every run, so re-running it
  IS the check-and-correct.
- **Resident/loaded check (consumer side, per frame):** the reader's
  `DftShieldingStore::hasJob(orig)` distinguishes "no job registered"
  from "registered but not yet parsed"; a frame that registers but fails
  validation is remembered absent and not retried
  (`model/DftShieldingStore.h:84-96`, `.cpp:46-71`). Rediscovery logs
  `DFT frames loaded=… gaps=…` (`rediscover/RunData.cpp:119,131`) and
  `dft_rows` (`rediscover/RunData.cpp:143-145`). Both surface a gap as a
  gap — never a faked value (memory `best_with_backup_and_locate`).

The CHECK is therefore "does the registered frame set equal the
completed-on-disk frame set, and does each registered frame validate."
The first is a set comparison; the second the consumers already do at
load. No new persistent index is introduced.

### MAKE UP THE DIFFERENCE (incremental, idempotent)

- **Mechanism:** re-run the registrar over the completed jobs
  (`lgs_write.py --force --dft-jobs-dir <campaign>/jobs <calcset_root>`).
  Because the emitted `dft.frames[]` is a pure function of the
  `orca_exit_code == 0` jobs on disk, re-running after the campaign grows
  produces a superset; nothing already-registered is lost; re-running
  with no new completions yields a byte-identical `.LGS`. That is the
  idempotency the brief asks for. Each subsequent reader-open / rediscover
  run then picks up the newly-listed frames for free, since both read
  `dft.frames[]` fresh at load.
- **Granularity:** registration is whole-frame. A frame is either listed
  (job complete + meta present) or not. There is no partial-frame state
  to reconcile; the loader's per-atom validation is all-or-nothing per
  frame. So "make up the difference" = "append the frames that finished
  since last registration," achieved by re-emit.
- **Multi-stride honoured:** the dense stride-1 extraction is unchanged
  across re-registrations; only `dft.frames[]` (the stride-2 list) grows.
  The original-index join keeps dense rows and sparse DFT aligned at
  every step.

Net: the flow is **(re)run `lgs_write.py` → both consumers re-read
`dft.frames[]` on next load**. The "difference" is made up at
registration time and consumed lazily; no consumer caches a stale frame
list across loads.

---

## 5. Reusable vs net-new

### Reusable as-is (the bulk — this is mostly a hygiene job)

- `tools/lgs_write.py` — the registrar: jobs-dir walk, typed
  `frame_index`, `orca_exit_code == 0` filter, sort, stride summary,
  `--force`, `--dft-jobs-dir` (`tools/lgs_write.py:145-219,334-454`).
- `CalcsetManifest::Load` + `DftFrame` — `.LGS` parse, `dft` block,
  per-frame path resolution + existence check
  (`io/CalcsetManifest.cpp:459-526`).
- `DftShieldingLoader::LoadAndValidate` + `OrcaShieldingParser` — the
  shared `.out` parse + validation, used identically by both consumers
  (`io/DftShieldingLoader.cpp`, `io/OrcaShieldingParser.cpp`).
- `DftShieldingStore` (reader) and `DftFrameSet` + `FrameMap` (rediscovery)
  — the two consumer-side homes for the parsed frames, both keyed by
  original index (`model/DftShieldingStore.*`, `rediscover/RunData.*`).
- The `orca_dft:*` catalogue entries and the strip-chart sampling path
  (`model/TrajectorySignalCatalog.cpp:1037-1061`).

### Genuinely new (small, and mostly operational/doc — no C++ required)

- **A canonical-`.LGS`-per-dir decision + cleanup.** Resolve Gap B: the
  calcset dir must hold exactly one `.LGS` (rename/retire the older
  `1p9j-calibration-with-dft.LGS`, or move the two datasets into separate
  roots), OR consumers must always be pointed at a specific `.LGS` file
  path, never the directory. This is a curation/ops choice, not code
  (the loader's one-`.LGS`-per-dir rule is correct and stays).
- **A re-registration trigger for the post-hoc cadence.** Today
  `lgs_write.py` is run by hand. The "incremental as the campaign grows"
  ask wants a defined moment to re-run it (e.g. a thin wrapper invoked
  before a rediscover run / before opening the reader on a live
  campaign, or a documented "after each consolidation, re-emit"). This is
  a one-line driver + a note, not new parsing. NOTE per
  `feedback_project_less_plastic_prefer_no_code`: prefer the documented
  manual/wrapper re-run over building a watcher; surface the choice to
  the lead before writing anything.
- **(Optional, only if the lead wants glance-level staleness in-app.)** A
  read-only coverage line in the reader's DFT log — "registered N of
  campaign_target M" — derived from values already on the manifest
  (`dft.frames.size()` vs `dft.campaign_target_frames`,
  already logged at `app/ReaderMainWindow.cpp:476-479`). Display polish,
  not mechanism; defer unless asked.

No schema change. No new map. No producer/extractor change. No H5 write.

---

## 6. Staged build plan (scope only)

Ordered, each stage independently verifiable, smallest-risk first. This
is a registration-hygiene + re-cadence task, not a code feature.

1. **Re-register the complete campaign (REGISTER).** Re-run
   `lgs_write.py --force --dft-jobs-dir /shared/2026Thesis/1p9j-orcas/jobs`
   against the calcset root so `dft.frames[]` covers all 751 completed
   even frames (0 → 1500), superseding the 660-frame (0 → 1318) list.
   Verify: emitted `frames[]` length == count of `orca_exit_code == 0`
   jobs; `frame_stride.last == 1500`. (Operational; no build.)

2. **Resolve directory ambiguity (Gap B).** Decide the one canonical
   `.LGS` for the calcset dir and retire/relocate the other so a
   directory-path load resolves uniquely (spec
   `CALCSET_MANIFEST.md:46-48`). Verify: `CalcsetManifest::Load(<dir>)`
   succeeds (no "multiple .LGS" error). (Curation; surface the
   keep/retire choice to the lead — `feedback_data_format_discuss_first`.)

3. **Confirm both consumers light up from the one registration (CHECK).**
   Reader: open the calcset; `orca_dft:{total,diamagnetic,paramagnetic}`
   sample non-empty across 0 → 1500 even frames; the DFT-store log shows
   `frames=751`. Rediscovery: `RunLoader::Load` logs `DFT frames
   loaded=751 gaps=0` and `dft_rows` == 751. Same `dft.frames[]`, both
   green. (Read-only verification; no new code.)

4. **Define the make-up cadence (MAKE UP THE DIFFERENCE).** Document (or,
   if the lead opts in, a thin wrapper) that re-runs `lgs_write.py
   --force` before a live-campaign rediscover run / reader open, so each
   load picks up newly-complete frames idempotently. Verify idempotency:
   re-run with no new completions ⇒ byte-identical `.LGS`. (Doc-first;
   code only on explicit opt-in.)

5. **(Optional) reader staleness line.** Only if requested: surface
   "registered N / campaign_target M" read-only in the DFT log. (Polish.)

The spine stays untouched: the typed model is populated from `.LGS` +
sidecar + `.out` exactly as today; the registration is just kept current
and unambiguous, and both consumers ride the one `dft.frames[]` map.

---

## Appendix — load-bearing citations

- Reader builds store iff `manifest.dft` present, from `dft.frames`:
  `app/ReaderMainWindow.cpp:467-469`.
- Reader DFT parse/validate (read-only, no H5):
  `io/DftShieldingLoader.cpp:45-111`; parser `io/OrcaShieldingParser.cpp:46-126`.
- Store is a release-after lazy source; history in strips not store / H5:
  `model/DftShieldingStore.h:18-38`, `.cpp:46-93`.
- T2 kept whole (raw 3×3 + T0/T1/T2; sample T0 or |T2|):
  `model/DftShielding.h:20-53`, `model/DftShieldingStore.cpp:73-86`.
- Catalogue advertises `orca_dft:*` unconditionally:
  `model/TrajectorySignalCatalog.cpp:250-251,1037-1061`.
- Rediscovery pulls the SAME `dft.frames[]` via the SAME loader:
  `rediscover/RunData.cpp:120-130`; join by original index +
  basis check `:46-74`.
- One-`.LGS`-per-dir rule (Gap B): `io/CalcsetManifest.cpp:166-177`;
  spec `CALCSET_MANIFEST.md:46-48,184-187`.
- Registrar (REGISTER + MAKE-UP): `tools/lgs_write.py:145-219,334-454`.
- `dft.frames[]` authoritative; `frame_stride` decoration; original-index
  join: spec `CALCSET_MANIFEST.md:142-178`.
- Campaign complete (751 even frames 0→1500): `/shared/2026Thesis/1p9j-orcas/jobs/`
  (751 dirs, last `…_f001500_t15000.0`), `status/jobs.csv` (header + 751).
- Dense LGS is populated, stale at 1318: 660 `frames[]`, even 0→1318,
  `frame_stride.last: 1318`, `campaign_target_frames: 751`
  (`…/1p9j-calibration-dense-mopac-live-orca.LGS`).
- Atom-count gate passes: meta `atom_count: 846`; `.out` has 846
  `Nucleus` in the `CHEMICAL SHIELDINGS` section; dense extraction
  topology atom == 846 (grounding doc §2).
