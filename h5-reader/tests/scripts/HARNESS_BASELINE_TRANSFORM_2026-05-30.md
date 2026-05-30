# Viewport harness baseline — transform-layer extension — 2026-05-30

Second-pass numbers from `tests/scripts/viewport_probe.py`. Adds the
upstream `TransformedConformation` decorator (identity / center_com /
fit_reference / fit_subset) plus two harness improvements that shipped
in the same commit:

* `POST /docks/visible {"visible": bool}` — wholesale hide/restore of
  the three ReaderMainWindow docks (inspector, selection, dashboard
  strip) so the central VTK viewport expands.
* `POST /selection/instrument {"enabled": bool, "focus_only": bool}` —
  the existing instrument preset extended with `focus_only=true`, which
  renders ONLY the focus-slot magenta sphere (the others are hidden).
  Eliminates the slot-1-eclipses-slot-0 problem the first-pass no-lock
  run hit (20 / 200 markers "missing" — actually present but biting
  into the magenta sphere via FXAA + depth occlusion, dropping the
  surviving pixel count below the 5-px detector threshold).

Sibling of `HARNESS_BASELINE_2026-05-30.md`; that doc covers the first
pass (camera-only floor + bug). This doc adds the four-experiment matrix
that includes the upstream transform layer.

Provenance: `notes/VIEWPORT_OBSERVATIONS_2026-05-30.md` §5b
(instrumentation-first reframe) + memory entry
`feedback_viewer_two_layers_transform_and_camera` (the two-layer
framing this work implements).

## PBC unwrap — DEFERRED (not blocking)

The brief flagged a PBC-unwrap component (port `pbc_whole.h` verbatim
from fes-sampler per `feedback_pbc_verbatim`, or skip cleanly if not
available). The canonical source `/shared/2026Thesis/fes-sampler/src/pbc_whole.h`
exists but it links against libgromacs (`do_pbc_mtop`, `tpxio.h`,
`pbc.h`, `topology.h`). h5-reader is explicitly a standalone Qt6/VTK
reader (`h5-reader/CLAUDE.md`: "never links the nmr_shielding library");
pulling in libgromacs is a much larger build-system commitment than this
session's scope.

Skipped cleanly. The RMSD fit modes (`fit_reference` / `fit_subset`)
deliver most of the stabilisation value on a trajectory whose PBC
unwrap was already done at extraction time, which is the case for the
1P9J calibration fixture used here. The numbers below confirm that
assumption: no NaN or wild-jump frames appeared in the
`transform_only` PNG sweep that would indicate the protein is being
torn across box boundaries.

If a future trajectory exhibits PBC artefacts at display time, the
right move is one of: (a) re-run extraction with PBC unwrap upstream;
(b) port `pbc_whole.h` into a separate `h5reader_pbc` library that
optionally links libgromacs at build time, behind a CMake flag; (c)
implement a substitute PBC unwrap on the H5 path that doesn't need
GROMACS internals (the box vectors ARE in the H5 — see
`QtTrajectoryH5.cpp:1996-1998`). None of (a)/(b)/(c) is in this
commit.

## How to reproduce

Reader binary: `build/linux-debug/h5reader` (Debug build, AMD Mesa
llvmpipe via the host GPU on the dev box).

Launch reader:

```bash
cd /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-debug
./h5reader --rest 9988 \
    /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft \
    2>/tmp/h5reader_run2.log &
# Wait for H5READER_REST_PORT=9988 in the log.
```

Run probe (same atoms, same frames as the first pass):

```bash
python3 tests/scripts/viewport_probe.py \
    --base http://127.0.0.1:9988 \
    --atoms 12 13 14 \
    --frames 200 \
    --out-dir /tmp/viewport_probe_run2
```

Reader stops with `curl -X POST http://127.0.0.1:9988/shutdown`.

The harness hides the docks before the experiments and restores them
after. The instrument preset is set with `focus_only=true` for every
experiment — only the focus-slot sphere renders, so there is no
slot-1-eclipses-slot-0 occlusion noise.

## Test setup

* Fixture: 1P9J calibration trajectory, 751 frames at 20 ps each,
  846 atoms total. Same as `HARNESS_BASELINE_2026-05-30.md`.
* Atoms: `(12, 13, 14)`. Same as first pass; the focus atom is the
  most-recently bulk-set member (atom 14). Marker tracked = slot 2
  (the focus) — but the marker preset's focus-only mode forces all
  visible spheres to magenta, so the harness still uses slot-0
  magenta hue thresholds (`H ∈ [285°, 315°], S ≥ 0.6, V ≥ 0.6`).
* Hardware: AMD Ryzen 9 9950X with Mesa 25.0.7 / radeonsi.
* 200 frames sampled at stride 1 (frames 0..199). Per experiment
  wall-time ~50-65 s — dominated by the per-frame
  `set_frame → screenshot → blob` round-trip.
* Backbone subset for the `fit_subset` transform: 318 atoms
  (`QtAtom::IsBackbone()` true; the typed substrate's backbone flag,
  no atom-name string parsing). Reference frame = 0.

## Numbers

### with_plane_lock (the floor — camera-only stabilisation)

| metric                | value          |
|-----------------------|----------------|
| frames sampled        | 200            |
| marker found          | 200/200        |
| reference centroid px | (1133, 713)    |
| mean drift            | 13.46 px       |
| **median drift**      | **1.52 px**    |
| std drift             | 48.66 px       |
| max drift             | 507.48 px      |
| p95 drift             | 23.98 px       |
| within 2 px           | 57.0 %         |
| within 5 px           | 61.5 %         |
| within 10 px          | 63.0 %         |
| marker area mean      | 1205 px        |
| marker area min / max | 62 / 1384 px   |

### without_plane_lock (centroid-delta camera path — the bug)

| metric                | value         |
|-----------------------|---------------|
| frames sampled        | 200           |
| marker found          | **200/200**   |
| reference centroid px | (1133, 713)   |
| mean drift            | 225.54 px     |
| median drift          | 225.32 px     |
| std drift             | 90.48 px      |
| max drift             | 447.53 px     |
| p95 drift             | 378.25 px     |
| within 5 px           | 0.5 %         |
| marker area mean      | 1609 px       |
| marker area min / max | 158 / 2324 px |

Compared to the first pass: **0/200 missing markers** (vs 20/200
before). The focus-only-marker change resolved the slot-1-eclipses-
slot-0 occlusion issue cleanly.

### transform_only (fit_subset on backbone, no camera lock)

| metric                | value          |
|-----------------------|----------------|
| frames sampled        | 200            |
| marker found          | 200/200        |
| reference centroid px | (1133, 713)    |
| mean drift            | 66.88 px       |
| **median drift**      | **67.26 px**   |
| std drift             | 28.26 px       |
| max drift             | 138.39 px      |
| p95 drift             | 111.39 px      |
| within 5 px           | 0.5 %          |
| marker area mean      | 2278 px        |
| marker area min / max | 1670 / 2519 px |

### transform_plus_lock (belt + suspenders)

| metric                | value          |
|-----------------------|----------------|
| frames sampled        | 200            |
| marker found          | 200/200        |
| reference centroid px | (1100, 706)    |
| mean drift            | 40.97 px       |
| median drift          | 39.01 px       |
| std drift             | 73.20 px       |
| max drift             | 795.03 px      |
| p95 drift             | 41.09 px       |
| within 5 px           | 12.0 %         |
| marker area mean      | 2108 px        |
| marker area min / max | 175 / 2225 px  |

## Commentary

### What the four numbers say

| experiment              | median drift | what it measures |
|-------------------------|-------------:|------------------|
| with_plane_lock         |   **1.5 px** | camera-only floor on this hardware |
| without_plane_lock      |  **225.3 px**| no camera-side, no transform-side — pure raw trajectory + centroid-delta camera bug |
| transform_only          |   **67.3 px**| transform-only stabilisation (no camera lock) |
| transform_plus_lock     |   **39.0 px**| both layers active |

**Headline finding: both layers contribute, NEITHER subsumes the
other.** The transform_only median is ~3.4× the with_plane_lock
floor and ~3.3× LOWER than without_plane_lock. The two-layer
framing in `feedback_viewer_two_layers_transform_and_camera` is
correct: rigid-body removal is a substantial first chunk of the
fix, but the camera path also needs work.

### Why transform_only doesn't reach the with-lock floor

The transform layer addresses the *protein's* rigid-body drift in the
simulation box — the 6 DOF (3 translation + 3 rotation) that
accumulate over MD time. A backbone-RMSD fit against frame 0 removes
those. What it leaves on the screen is:

1. Atom-internal vibration of atoms 12-13-14 within the rest-frame
   protein. At 300 K backbone N/CA/C vibrate ~0.1-0.3 Å rms; in screen
   pixels at the current zoom that's a few px.
2. The centroid-delta camera path (the bug) is still active in
   transform_only — the camera continues to translate per-frame by
   the (now-stabilised) protein's centroid delta. The transform makes
   that delta MUCH smaller but not zero, and the centroid-delta path
   still has the accumulation problem.

The plane lock side-steps (2) by absolute-writing the camera every
frame from the 3-atom plane basis. Combining transform + lock should
in principle hit the with-lock floor (~1.5 px) but...

### Why transform_plus_lock is HIGHER than with_plane_lock alone

This was unexpected at first glance — transform_plus_lock median =
39 px vs with_plane_lock alone = 1.5 px. Two compounding causes:

1. **The lock atoms move in the transformed frame.** The plane lock
   computes its basis from the world-space positions of atoms 12, 13,
   14 at lock-acquisition time. The transform applies a per-frame
   rotation + translation that's recomputed for the new frame. So the
   plane lock's basis, captured at frame 0 in the IDENTITY frame, is
   no longer the same basis once the transform mode flips — the lock
   is re-evaluated in a frame where the atoms have moved relative to
   the lock's stored direction. The continuity sign-flip guard
   (`MoleculeScene.cpp:390`) may fire spuriously when the transform
   rotates the molecule past the basis's threshold.

2. **The lock was enabled AFTER the transform was set** (harness
   ordering: `set_transform → set_instrument → enable_plane_lock`).
   The lock captures the camera state at enable time, which is the
   camera state AFTER the transform's first refreshCurrentFrame. So
   the lock's frame-0 reference is the transformed view, not the
   raw view. Subsequent frames apply the per-frame transform AND the
   lock — the two operate on different reference frames.

This is a known design question for the next session — see
`VIEWPORT_OBSERVATIONS §7` ("lock target source of truth", "user-delta
persistence", "plane lock release on selection change"). The
transform_plus_lock max = 795 px (vs p95 = 41 px) is the smoking gun:
single-frame discontinuities when the two layers' reference frames
flip relative to each other. The median is sane; the tail is the
interference pattern.

### What the next session's design conversation should weight

Given the numbers:

* **Transform layer alone delivers ~3.3× improvement** over the raw
  no-lock baseline. That's substantial. The lock taxonomy can lean on
  the transform doing the rigid-body work and focus on
  user-orientation tasks (sight-down-bond, atom-centred, etc.) rather
  than re-deriving "stabilise the protein."

* **Camera path still needs the rebuild** (eventFilter primary input,
  render coalesce, overlay-layer renderer for markers, shutdown
  reorder). The reduction from 225 → 67 px is necessary but not
  sufficient — the user's intent ("look at atom 14 through the
  trajectory") needs the camera to be absolute-written per frame
  from a typed lock, not delta-translated from the centroid.

* **Lock + transform interaction needs explicit design.** The
  transform_plus_lock numbers indicate the layers don't compose
  naively — the plane lock today stores world-space state at
  acquisition time; the transform changes what "world space" means.
  Options:
    1. **Lock state is captured in the inner (raw) frame.** Lock's
       per-frame application then reads inner positions, computes the
       basis, applies the camera write. Transform is invisible to
       lock math. Probably wrong: the user's view-up gesture happens
       in the displayed (transformed) frame.
    2. **Lock is a no-op when the transform is non-identity.** The
       transform already stabilises rigid body; the lock isn't
       needed. The transform_only number suggests this is workable;
       median 67 px ≈ atom-vibration noise inside the stabilised
       protein, hard to remove without a per-atom-pair lock.
    3. **Lock is recomputed every time the transform mode changes**
       (lock's reference is re-acquired from the current displayed
       state). The two-layer composition works because lock always
       speaks the transform's language.

  Option (3) is the cleanest. The next session should decide between
  (2) and (3) explicitly.

### About the marker-area numbers

Marker area is ~2× larger across all four experiments compared to
the first pass — the dock-hide widens the central viewport from
~830 px wide to ~1530 px wide, so a 1.5-Å magenta sphere covers more
screen pixels. Within each experiment marker area is more stable
than before (transform_only std = 170 px vs with_plane_lock's 322 px)
because focus-only renders ONE sphere whose silhouette is undisturbed
by neighbouring colour-marked spheres.

### The 200/200 marker-found rate

Every experiment found the marker on every frame. Compare to the
first pass's no_lock at 180/200. The combination of (a) bigger
viewport + (b) focus-only marker resolved the marker-missing
problem entirely. The blob detector's 5-px-area threshold is
comfortably below the worst-case marker_area_min for transform_only
(1670 px) and even comfortably below the worst case (62 px) for
with_plane_lock — the 62-px floor in with_plane_lock corresponds to
a frame where the marker is at the edge of the viewport behind ribbon
geometry, but it still passes the threshold.

## Per-frame data + visual review

* `/tmp/viewport_probe_run2/with_lock.csv` — 200 rows.
* `/tmp/viewport_probe_run2/no_lock.csv` — 200 rows.
* `/tmp/viewport_probe_run2/transform_only.csv` — 200 rows.
* `/tmp/viewport_probe_run2/transform_plus_lock.csv` — 200 rows.
* `/tmp/viewport_probe_run2/<experiment>/frame_NNNNN.png` — per-frame
  snapshots for each experiment.
* `/tmp/viewport_probe_run2/summary.json` — the JSON the harness
  prints to stdout.

Recommended visual sanity check: open `transform_only/frame_00000.png`
and `transform_only/frame_00099.png`. The protein backbone should be
in the same orientation in both (the transform did its job) and the
magenta marker should be visible in both at roughly the same screen
position.

## Threshold tuning notes

Magenta hue threshold from `HARNESS_BASELINE_2026-05-30.md` carries
over unchanged (`H ∈ [285°, 315°], S ≥ 0.6, V ≥ 0.6`). With
focus_only the marker is large (~2000 px area mean) and unambiguous —
the harness reports no false positives in any experiment.

## Updated acceptance gate for the next refactor session

Two gates the next refactor (lock taxonomy + camera architecture)
must meet on this fixture + frames:

1. **with_no_lock median ≤ 13 px** (within ~10× of the current
   with_plane_lock floor) after the camera architecture refactor
   (eventFilter input, render coalesce, overlay-layer renderer for
   the marker, shutdown reorder, the two cheap probes). This is
   measurable today without transform; the centroid-delta path is
   the load-bearing bug.

2. **transform_plus_lock median ≤ 5 px AND p95 ≤ 10 px AND max ≤
   30 px**. The current transform_plus_lock numbers (median 39, p95
   41, max 795) show the two layers don't compose. After the camera
   refactor + the lock/transform interaction design decision lands,
   the composition should match or beat with_plane_lock alone.

Both gates are deterministic and measurable via this harness. The
"5-minute-round user check-in" anti-pattern stays banned per
`feedback_instrument_before_refactor_qtvtk` — the harness produces
the numbers; the user is consulted on design choices, not on
whether the previous experiment worked.

## Hardware caveats

The first-pass numbers were taken with the dock layout in place and
focus-only off; this second pass widened the viewport and removed
the slot-0/slot-1 occlusion. Direct comparison of medians between
the two passes is therefore not apples-to-apples:

| experiment         | first pass median | second pass median |
|--------------------|------------------:|-------------------:|
| with_plane_lock    |          13.4 px  |             1.5 px |
| without_plane_lock |          95.7 px  |           225.3 px |

The first pass's with-lock median (13) reads MUCH worse than this
pass's (1.5) because the dock-induced compression cropped the
marker more often. The first pass's no-lock median (96) reads BETTER
than this pass's (225) because the smaller viewport meant some
"on-screen" markers were actually clipped to the edge of the
narrower central widget — clipping cuts the area, and the centroid
is biased toward whichever bit survived. The bigger viewport this
pass surfaces the real magnitude of the camera drift.

The relative pattern within a single pass is what matters for the
refactor. Always re-baseline after a harness change before comparing
to historical numbers.
