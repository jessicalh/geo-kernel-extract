# Viewport harness baseline — 2026-05-30

First-run numbers from `tests/scripts/viewport_probe.py` against the
1P9J calibration fixture. Captures the "marker drift floor" (camera
plane lock owns the camera) and the "floor + current bug" (centroid-
delta camera-follow path owns the camera). The next session's viewport
refactor (eventFilter primary input, render coalesce, lock taxonomy,
shutdown reorder) targets the no-lock numbers — without breaking the
with-lock numbers.

> **Supersession**: a second pass added the upstream
> `TransformedConformation` layer (identity / center_com /
> fit_reference / fit_subset), the `POST /docks/visible` endpoint, and
> the single-focus-atom marker mode. See
> [`HARNESS_BASELINE_TRANSFORM_2026-05-30.md`](HARNESS_BASELINE_TRANSFORM_2026-05-30.md)
> for the four-experiment matrix on the same fixture + atoms. That
> doc updates the acceptance gate for the lock-taxonomy + camera
> refactor session.

Provenance for these numbers: `notes/VIEWPORT_OBSERVATIONS_2026-05-30.md`
§5b (instrumentation-first reframe), §4 (paint-some-frames probe), §1
(why it bounces).

## How to reproduce

Reader binary: `build/linux-debug/h5reader` (Debug build, AMD Mesa
llvmpipe via the host GPU on the dev box).

Launch reader:

```bash
cd /shared/2026Thesis/nmr-shielding/h5-reader/build/linux-debug
./h5reader --rest 9988 \
    /shared/2026Thesis/shielding-calcsets/data/trajectories/1p9j-calibration-with-dft &
# Wait ~5s for Qt + VTK + topology load. Check stderr for
# H5READER_REST_PORT=9988 then GL: vendor / renderer / version.
```

Run probe:

```bash
python3 tests/scripts/viewport_probe.py \
    --base http://127.0.0.1:9988 \
    --atoms 12 13 14 \
    --frames 200 \
    --frame-stride 1 \
    --out-dir /tmp/viewport_probe_run1
```

Reader stops with `curl -X POST http://127.0.0.1:9988/shutdown`.

## Test setup

* Fixture: 1P9J calibration trajectory, 751 frames at 20 ps each,
  846 atoms total.
* Atom triple `(12, 13, 14)`, hard-coded as the first stable backbone
  triple inside the protein (positions at frame 0 in Å):

  | atom | x       | y       | z       |
  |------|---------|---------|---------|
  | 12   | 88.557  | 67.853  | 20.056  |
  | 13   | 88.132  | 68.821  | 19.791  |
  | 14   | 87.770  | 67.280  | 20.547  |

  Pairwise distances: d(12,13) = 1.09 Å, d(13,14) = 1.75 Å,
  d(12,14) = 1.09 Å. Plane normal magnitude ≈ 1.14 (non-degenerate;
  `/plane-lock/enable` accepts the triple). Geometric meaning: bond
  to non-bond pair — a stable plane that survives MD vibration.
* Marker tracked: slot 0 (pure magenta `#FF00FF`). Hue threshold
  `H ∈ [285°, 315°], S ≥ 0.6, V ≥ 0.6` against the
  RGB-→-HSV transform of the PNG.
* Hardware: AMD Ryzen 9 9950X with on-die Mesa llvmpipe rasterising
  via the host's iGPU path (`OpenGL renderer = AMD Ryzen 9 9950X
  16-Core Processor (radeonsi, raphael_mendocino, …)`,
  Mesa 25.0.7). FXAA on, depth peeling off (per
  `MoleculeScene.cpp:80-81`).
* 200 frames sampled at stride 1 (frames 0..199 of the 751 available).
  Per experiment: ~30 s wall time, dominated by the per-frame
  `set_frame → screenshot → blob` round-trip.

## Numbers

### With plane lock (the floor — what we want to preserve)

| metric                | value      |
|-----------------------|------------|
| frames sampled        | 200        |
| marker found          | 200/200    |
| reference centroid px | (748, 415) |
| mean drift            | 45.27 px   |
| median drift          | 13.38 px   |
| std drift             | 93.01 px   |
| max drift             | 407.03 px  |
| p95 drift             | 290.47 px  |
| within 2 px           | 11.0 %     |
| within 5 px           | 28.0 %     |
| within 10 px          | 41.5 %     |
| marker area mean      | 152.8 px   |
| marker area min / max | 6 / 303 px |

### Without plane lock (centroid-delta path — the bug)

| metric                | value          |
|-----------------------|----------------|
| frames sampled        | 200            |
| marker found          | 180/200        |
| **marker MISSING**    | **20/200**     |
| reference centroid px | (721, 413)     |
| mean drift            | 111.81 px      |
| median drift          | 95.65 px       |
| std drift             | 80.46 px       |
| max drift             | 338.46 px      |
| p95 drift             | 317.40 px      |
| within 2 px           | 0.6 %          |
| within 5 px           | 0.6 %          |
| within 10 px          | 1.1 %          |
| marker area mean      | 140.3 px       |
| marker area min / max | 5 / 377 px     |

## Commentary

The bug is measurable. The signal-to-noise is roughly 7× on the median
and 2.5× on the mean — comfortably larger than per-frame jitter and
larger than the run-to-run reproducibility noise that GPU drivers can
inject. The next session can demand "median drift under no_lock comes
down to within 2× of the with_lock median" as the acceptance gate for
the refactor and have a deterministic test for it.

Three observations that the next session's refactor needs to land on:

* **20 frames out of 200 have the marker entirely off-screen** in
  no_lock mode. The atom hasn't moved that far in real space — the
  protein is at ~300 K vibrating a few tenths of an Å — so the camera
  has drifted away from where the atom is. This is the centroid-delta
  accumulation bug described in `VIEWPORT_OBSERVATIONS §1.1`:
  per-frame deltas keep adding (or worse, PBC wraparound + a single
  teleport delta) until the atom leaves the frustum.

* **With plane lock, median = 13 px not 0.** The plane lock pins the
  camera's orientation to the basis of the (12, 13, 14) plane, but
  atom 12 itself moves within that plane as the molecule vibrates;
  the marker centroid moves with it. 13 px ≈ 0.4 Å in the rendered
  scene's pixel scale, which is consistent with atom-scale thermal
  motion. The plane-lock max (407 px) reflects the normal-sign
  continuity flip handled by `applyCameraPlaneLock` line 390 —
  when the protein crosses the configuration that flips the
  cross-product sign, the camera rotates 180° around the plane's
  origin and the marker lands on the other side of the screen. The
  algorithm is correct; the marker drift reads it as a transient
  "drift" pulse. The median is the right central tendency, not the
  mean.

* **Marker area variance is high in both modes** (min 5–6 px, max
  300–377 px). The marker is occluded by molecule geometry on some
  frames (the imposter spheres of the atoms in front bite into the
  marker's silhouette). The blob finder still picks the largest
  connected component, but its centroid measurement is biased toward
  whichever bit of the marker survived the occlusion. For the lock-
  taxonomy work this is fine — we're measuring the relative drift
  between two modes on the same atom — but it is the floor that the
  overlay-layer-renderer refactor (instrument-mode marker in its own
  `vtkRenderer::SetLayer(1)`, depth-occlusion-immune) would push down
  toward the GPU rasteriser noise level (~1 px).

## Per-frame data + visual review

* `/tmp/viewport_probe_run1/with_lock.csv` — 200 rows: frame,
  centroid_x, centroid_y, area_px, png_sha1.
* `/tmp/viewport_probe_run1/no_lock.csv` — same, with blank rows where
  the marker was missing.
* `/tmp/viewport_probe_run1/with_lock/frame_00000.png` ..
  `frame_00199.png` — per-frame snapshots.
* `/tmp/viewport_probe_run1/no_lock/frame_00000.png` .. ditto.
* `/tmp/viewport_probe_run1/summary.json` — same JSON the script
  prints to stdout.

Recommended visual sanity check: open `frame_00000.png` and
`frame_00100.png` for both experiments. With-lock the marker is at
roughly the same screen position; no-lock the marker is somewhere
else entirely (or invisible because it left the viewport).

## Threshold tuning notes

The magenta hue threshold (`H ∈ [285°, 315°], S ≥ 0.6, V ≥ 0.6`) was
chosen by visual inspection of the first snapshot under instrument
mode and a small loop checking pixel counts per slot. On the Mesa
llvmpipe path with FXAA on, the rendered marker stays well within
this range; ~250 pixels in the centre, dropping to ~50 when partially
occluded. No threshold tuning was required after the first run.

On other GPU paths the hue may shift slightly (driver-specific FXAA
edge desaturation, OpenGL gamma handling); re-tune by adjusting
`MARKER_HUES[0]` in `viewport_probe.py` if the marker hue lands
outside the band on a different host.
