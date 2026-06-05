# Basic-screen state + plan — 2026-06-05

## Where we are
- The reader's selected-metrics coherence + the "strips show one dot" bug are DONE and committed (`141e7b6`):
  single-source-of-truth selection controller, availability/renderability model roles, one capability table,
  selector reflects the picked row, dock reveal-on-add/run-load-reset, strip-history persistence across
  rebuilds, and the NPY snapshot-load fix (`TransformedConformation::loadSnapshot` now triggers the inner load
  instead of reading a stale resident slot). Verified value-exact vs independently-read NPY ground truth, both
  H5 and NPY backends. Fixed build is up on `:1` for the lead to mouse-confirm.
- **Pivot (lead, 2026-06-05):** stop adding; get the BASICS consistently right. Strip the fancy.

## The goal — a basic, consistent screen
One screen where you can:
1. **Pick one of two stabilisation modes** (and clearly see which you're in / switch between them).
2. **Zoom / pan / maneuver and KEEP that view** — the camera holds its spot; it does not reset or get hijacked.
3. **Select an atom** reliably.
Nothing else. No residue-only screen, no multi-select camera-lock modes, no half-working display options.

## The two stabilisation modes (kept) — what they mean
1. **Locked backbone** — the INDUSTRY STANDARD: least-squares (Kabsch) fit of each frame's **backbone** atoms
   onto a reference frame, removing global translation + rotation. The backbone holds visually still while
   sidechains / loops show their real motion. Cf. GROMACS `gmx trjconv -fit rot+trans` (Backbone group), VMD
   RMSD Trajectory Tool, PyMOL `intra_fit`, MDAnalysis fit transforms.
2. **Kabsch-stabilised "with give"** — removes tumbling/rotation but lets ALL real internal motion show (the
   backbone still flexes; nothing is hard-held). The misbehaving variant the lead calls the **"Kabsch
   gyroscope"**: a residual spin/precession, most likely because an all-atom fit of a flexing molecule has no
   stable frame. Investigate whether anchoring the fit to a stable subset removes the precession while keeping
   "give."
- **The switch is currently NOWHERE CLEAR** — there is no visible indicator of which mode is active or how to
  toggle. The basic screen must make this obvious.

## Temporal stabilisation — THE key refinement (lead decision, 2026-06-05)
The "Kabsch gyroscope" is NOT fixed by changing the fit target. Keep BOTH modes and add **temporal
stabilisation** of the alignment — the lead's instinct, and it's the right one:
- The alignment transform (the R/T mapping each frame onto the reference) is computed **per-frame,
  independently** — zero temporal coherence (inventory-confirmed: each frame Kabsch-fits independently to the
  cached reference, `TransformedConformation.cpp:157`). That per-frame jitter/wander IS the bad gyroscope; the
  actual tumbling is already removed by the fit.
- **Fix:** smooth the transform *sequence* over time — quaternion-smoothed rotation + smoothed translation over
  a rolling window — then apply the smoothed transform to the raw atoms. Orientation moves smoothly (rotation +
  tumbling-seen-over-time removed) while every atom's raw internal motion is preserved (the "give").
- **Mode-agnostic:** both locked-backbone and with-give get the smoothing; the ONLY difference between them
  stays the fit target (backbone subset vs all atoms).
- **One knob:** smoothing window length / time constant — sane default, tunable, **dialled in visually on `:1`**
  (too little → still jitters; too much → slow real motion lags).
- **Implementation shape:** precompute per-frame Kabsch (subset → reference frame 0) across the loaded
  trajectory, temporally smooth the (R,T) sequence, `atomPosition(f) = R'_f · raw_f + T'_f`. Replaces the
  per-frame-lazy cache in `TransformedConformation`.

## Inventory outcome (codex `bulfdtt7m`, READ-ONLY) → `notes/DISPLAY_MODE_INVENTORY_2026-06-05.md`
Cited `file:line` + VTK class-doc URLs. Confirmed mappings: **locked backbone = `FitSubset(BackboneSubset)`**,
**with-give = `FitReference`** (all-atom). **Free camera is sound and HOLDS the view** (no-op per-frame write).
Atom select model + marker are sound + index-stable, but the **pick is by ray-distance not frontmost/depth**
(can pick a hidden atom; also truncates click coords) — a real reliability gap. Only 4 static modes render;
~14 `static.*` are advertised-but-dead; `strip.spectrum`/`static.tensor` half-wired. All judgments accepted.

## Choices taken (lead, 2026-06-05)
- **Exactly two stabilisation modes**: locked backbone + Kabsch-with-give. No others.
- **No multi-select camera-lock modes** (Focus / Newman / Plane lock / Atom / Bond / Subset) — once the two
  stabilisation modes + a free camera that holds its view are in place, these aren't needed.
- **No residue-only view, no "fancy stuff."** Basics, consistent.
- **Free camera must persist** the user's zoom/pan/orbit; selecting an atom must work consistently.
- **Understand before cutting** — codex inventories everything; we weigh its judgments, don't rubber-stamp.

## Build plan (approved 2026-06-05) — in order
0. ✅ **Inventory done** (above). All judgments accepted; lead refined the gyroscope fix → temporal stabilisation.
1. **Temporal stabilisation** of the alignment (the spine) — both modes. Tune the smoothing window on `:1`.
2. **Named two-mode switch** ("Locked backbone" / "Kabsch with give") + a visible current-mode label
   (replaces the lone unlabeled "All-atom fit" checkbox, `ReaderMainWindow.cpp:887`).
3. **Reveal = highlight only** (no camera-mode hijack, `MoleculeScene.cpp:549`) + **depth-aware atom pick**
   (frontmost among hits, no int-truncation of click coords, `QtAtomPicker.cpp:55`).
4. **Cut** the camera locks (Focus/Newman/Plane/Atom/Bond/Subset), the Metrics/dashboard surface + the ~14
   advertised-but-dead `static.*` modes + half-wired `strip.spectrum`/`static.tensor`, the Butterfly/B-field
   overlays, the Ribbon/Rings toggles, the hidden harness marker, the dead `SignalBuffer::isValidAt`. (Keep
   the dashboard behind an Advanced path only if needed later — not deleted, just off the basic screen.)

Each step: codex grinds (qt skill loaded) → I verify on `:1` (off the lead's hands) → lead mouse-confirms.
Floor before this work: `141e7b6`.

## Resolved (were open questions)
- locked backbone = `FitSubset(BackboneSubset)`; with-give = `FitReference` (all-atom). [confirmed in inventory]
- The gyroscope = per-frame alignment jitter (all-atom-fit bias amplifies it); no temporal feedback in the
  code. FIX = temporal stabilisation (above), NOT a fit-target change. The reveal-click camera hijack (D1) is a
  SEPARATE issue, also fixed in step 3.
- Free view is held by `writeFree()` being a no-op; what overrides it = the camera LOCK modes + reveal-click,
  both cut/fixed in steps 3-4.

## Pointers
- Stabilisation/transform: `src/model/TransformedConformation.cpp`, `src/app/FitTargetMath.h`, transform
  actions in `ReaderMainWindow.cpp`. Camera: `CameraComposer`, `CameraInputFilter`. Selection: `QtAtomPicker`,
  `AtomSelection`. Display modes: `DisplayModeCapability.h` + `DashboardDisplayController` + the `*Panel`s.
- Punch list (broader): `notes/DESK_READY_PUNCHLIST_2026-06-05.md` (this work = its "Week sequencing" step 2).
- 2026-06-05: landed temporal smoothing for both display stabilisation modes by precomputing per-frame Kabsch transforms and smoothing only `(R,T)`; default half-window is 8 frames (17-frame window), `window=0` restores raw per-frame fits, and headless tuning is `POST /transform/smoothing {"window": int}` with `GET /transform` reporting `window`.

## SESSION END STATE (2026-06-05, ~98% ctx) — resume EXACTLY here
- **Committed `141e7b6`:** the selected-metrics coherence + the "strips show one dot" fix (top of doc).
- **UNCOMMITTED in the working tree (on top of `141e7b6`):** the **temporal stabilisation** above
  (`TransformedConformation` precompute + quaternion-smoothed `(R,T)`, REST `/transform/smoothing`). Build green.
  NOT committed — successor should review + commit with the lead's nod.
- **TUNING FINDING (lead, on `:1`, all-atom mode):** at **heavy smoothing (window=32) the molecule "duck walks"
  across the screen** — the windowed smoothing LAGS the true motion → a drift/waddle. That is the over-smoothing
  failure mode, and it is the argument for the next step (windowed smoothing of frame-0 fits inherently lags;
  the iterative-mean reference does not).
- **THE next step — iterative-mean reference (lead's "find the right space", APPROVED, was held until she'd
  tried the windowed version — she has):** we have the WHOLE trajectory in memory, so don't fit to an arbitrary
  frame 0 + smooth. Compute the optimal reference by **iterative superposition to the mean structure** (align
  all → mean → re-align → re-mean until converged), fit each frame to that converged mean. Global solution (no
  windowing → no duck-walk lag), unbiased (no frame-0 bias). Fit on a **stable subset** (backbone / rigid core)
  to the mean to kill the gyroscope contamination at the source; keep the smoothing window as an optional
  jitter knob. This is the truthful "remove tumbling & rotation ONLY": a uniform rigid transform cannot distort
  internal motion (preserves all distances/angles) — the only question is the orientation removed, and a
  stable-subset fit to the mean removes only the global rigid-body motion, so a real internal/domain rotation
  SURVIVES (that survival is the acceptance test).
- **Then:** (2) named two-mode switch ("Locked backbone" / "Kabsch with give") + visible current-mode label;
  (3) reveal = highlight-only (no camera hijack) + depth-aware atom pick (frontmost, no coord truncation);
  (4) cut the camera locks / Metrics surface / ~14 dead `static.*` modes / Butterfly+B-field / Ribbon+Rings
  toggles / hidden harness marker / dead `SignalBuffer::isValidAt`. (See DISPLAY_MODE_INVENTORY keep/cut/fix.)
- **WORKFLOW GAP — screenshots:** stabilisation tuning must be judged visually. The reader has `/screenshot`
  (`target` = `scene` / `window` / `widget`+`object_name`). Lead drives `:1` by mouse; the successor should run
  a REST instance on isolated **Xvfb `:99`** (software GL works: `LIBGL_ALWAYS_SOFTWARE=1 GALLIUM_DRIVER=llvmpipe`)
  and screenshot the scene across frames to assess stabilisation headlessly, and/or screenshot the lead's `:1`
  instance. This session tuned blind (no scene screenshots) — fix that first next time.
- **Codex briefs this session** (every fire prepends `notes/CODEX_PREAMBLE_QT_VTK.md` — standard qt-skill +
  VTK-docs-online-license preamble; KEEP doing this): INSTRUMENT_DASHBOARD, SELECTOR_REFLECTION,
  CAPABILITY_TABLE, SIGNALMODEL_ROLES, SELECTION_CONTROLLER, DOCK_LIFECYCLE, STRIP_SERIES_DUMP,
  STRIP_HISTORY_PERSIST, REReader_INVENTORY, TEMPORAL_STABILISATION. The inventory return is
  `notes/DISPLAY_MODE_INVENTORY_2026-06-05.md` (cited `file:line` + VTK doc URLs) — **READ IT, don't trust the
  summary.**
- **Known small problems:** `gromacs_energy` NPY cols 43 vs catalog 42 ("trusting the NPY", benign warning);
  dead `SignalBuffer::isValidAt`; `static.spectrum.power` capability asymmetry (renderable but picker-hidden);
  stale `static.tensor` comments claiming the glyph fires.
- **Live now:** reader on `:1` REST port in `/tmp/p1`, all-atom mode, window=32 (duck-walking). `141e7b6` floor.

## DUCK-WALK DIAGNOSED + TWO-PART FIX APPROVED (2026-06-05, new session)
**Diagnosis (measured on `:99`, 1P9J `1p9j-calibration-with-dft.LGS`, all-atom fit):** the duck-walk
is a TRANSLATION artifact, not a rotation one. Centroid of the rigid body across the run:
- `window=0`: centroid pinned at ≈(74.5, 48.4, 28.2), `|c|`≈93 Å, wanders <1 Å over all 751 frames.
- `window=16`: centroid SWINGS — by frame 750 it's at (94.0, 48.7, **−8.3**), `|c|`=106 Å (~40 Å excursion;
  molecule radius ~18 Å → walks off-frame). Screenshots `/tmp/ducks/w16_f750.png` (clipped at right edge)
  vs `/tmp/ducks/w0_f750.png` (centred) confirm visually.
**Root cause:** the raw GROMACS coords sit `|c|`≈93 Å from the origin and the transform rotates ABOUT THE
ORIGIN (`R·raw + T`). At window=0, R and T are a matched Kabsch pair → centroid lands exactly on the
reference. Smoothing R and T **independently** breaks the pairing; a small mismatch × the 93 Å lever arm =
a huge excursion. Worst at the trajectory ends, where the symmetric window clamps one-sided.
**Second, quieter instability:** even window=0 backbone isn't truly stable because frame 0 is an arbitrary
reference snapshot.

**Approved fix (lead, 2026-06-05) — two parts:**
- **A. Centroid-pin (kills the walk).** Apply `R'·(raw − cc_f) + cr` where `cc_f` = current-frame fit-set
  centroid, `cr` = reference centroid. Equivalently: smooth ONLY the rotation, then DERIVE the translation
  `T = cr − R'·cc_f` (never smooth T independently). Centroid maps to `cr` exactly every frame; lever arm
  drops 93 Å → ~18 Å. This is just the correct form — right regardless.
- **B. Iterative converged-mean reference (kills the jitter) — the lead's "use the full trajectory."**
  Replace the frame-0 reference with the iterative mean structure over the whole run (align all → average →
  re-align → converge), fit on the stable backbone core. Representative reference → far less per-frame
  jitter → little/no smoothing needed → no lag, no walk.
- Net: centroid pinned exactly, orientation from a stable-subset fit to the converged mean, ONE rigid
  transform per frame (internal motion preserved → a real domain rotation survives = the truth test).
  Smoothing window becomes an optional light knob, default likely **0**.

**In flight:** prior-art research codex (`btqus20it`, `notes/CODEX_BRIEF_STABILISATION_PRIORART_2026-06-05.md`
→ writes `notes/STABILISATION_PRIORART_2026-06-05.md`) — what PyMOL/VMD/GROMACS/MDAnalysis/Theseus do
(iterative-mean reference, centroid-centering, ML/stable-core weighting for the gyroscope, whether anyone
temporally smooths). Findings fold into the build brief before firing the build codex on `TransformedConformation`.
**`:99` reader live:** `build/linux-rwdi/h5reader`, REST port **36669** (software GL on Xvfb `:99`). Stay off `:1`.

## BUILT + VERIFIED (2026-06-05) — stabilisation fix landed (UNCOMMITTED)
Codex `blql3k3uu` rewrote `TransformedConformation.{h,cpp}` per `notes/CODEX_BRIEF_STABILISATION_BUILD_2026-06-05.md`
(grounded in `notes/STABILISATION_PRIORART_2026-06-05.md`): iterative converged-mean reference
(MDAnalysis `iterative_average`: max_iter=20, eps=1e-4 Å), centroid-pinned smoothing (smooth R only,
derive `T = cr − R·cc_f`), default `kDefaultStabilisationWindow = 0`. Build green (`ninja: no work to do`).
Also touched: `RestServer.{cpp,h}` (comment/doc updates only), `ReaderMainWindow.cpp` (one startup comment).
**`equiv_t2_e3nn.py` in the working tree is PRE-EXISTING rediscover work, NOT this change — leave it.**

**Opus adversarial review (`a92dd5d9`): SHIP.** No CRITICAL/MAJOR. Centroid invariant + truth bar both hold;
independent-T averaging gone; `cc_f` is the current-frame centroid; `cr` constant; smoothed-R/derived-T matched.
MINOR/NIT only: (a) no non-convergence warning in the iterative mean; (b) `rebuildTransformCache`/
`setStabilisationWindow` heavy path un-timed (brief asked for a timing log); (c) raw fits recomputed on every
window change (perf NIT); (d) stale line-range comment `FitTargetMath.h:13`. None blocks; fold (a)/(b)/(d) into
the pre-commit cleanup.

**Empirical proof on `:99` (port 33287, 1P9J):**
- EXACT centroid invariant — backbone fit-set (318 atoms) centroid is **constant to ~1e-13** across
  frame ∈ {0,375,750} and window ∈ {0,16}: `[73.84739424447596, 48.52753477276494, 28.39069193264223]`
  every time. The duck-walk is now mathematically impossible, not tuned away.
- All-atom window=16: the 9-atom sample centroid the OLD build swung to (94, 49, −8), `|c|`=106 now stays at
  (74–75, 48, 28), `|c|`=93–94 across the whole run (<1.5 Å real residual, no 40 Å artifact).
- Iterative mean converges in **4 iterations** (backbone 15 ms / all-atom 846-atom 36 ms) — no GUI stall.
- Screenshots `/tmp/ducks/new_w16_f{0,750}.png` + `new_bb_w0_f{0,750}.png`: molecule stays centred/framed
  (vs OLD `/tmp/ducks/w16_f750.png`, clipped off the right edge).

**To mouse-confirm:** the lead restarts her `:1` reader to pick up `build/linux-rwdi/h5reader` (the new binary).
The default opens as backbone fit, window 0. UNCOMMITTED — the lead owns git.
- 2026-06-05: landed TransformedConformation iterative-mean fit reference plus centroid-pinned smoothing; default stabilisation window is now 0. Convergence params: max_iter=20, eps=1e-4 A, seeded/anchored by `reference_frame` (default 0). REST shape unchanged (`POST /transform/smoothing {"window": int}` kept as optional rotation-only smoothing); UI unchanged.
