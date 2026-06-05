# Codex — TransformedConformation: iterative-mean reference + centroid-pinned smoothing

## (Standard qt6-cpp + VTK-docs preamble is prepended above — honour it. Verify from the code before changing anything; cite file:line.)

## Stakes / why this matters
The reader's whole job is to show *internal* molecular motion with the box tumbling/drift
removed. The current stabiliser has a measured bug: with temporal smoothing on, the whole
molecule walks ~40 Å across the screen and off-frame ("duck walk"). We diagnosed it on a live
trajectory and researched what mature tools do (PyMOL/VMD/GROMACS/MDAnalysis/Theseus). This
brief turns that into the fix. Get it right and the lead can finally hand the reader to her
advisor. Read `notes/STABILISATION_PRIORART_2026-06-05.md` and
`notes/BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md` (the "DUCK-WALK DIAGNOSED" section) first.

## The diagnosis, precisely (verify it yourself)
- Raw GROMACS coords sit ~93 Å from the origin. The per-frame Kabsch in
  `ComputeSubsetTransform` (`src/app/FitTargetMath.h:189-228`) is **already centroid-correct**:
  it demeans both point sets and sets `T = cr - R*cc`. That is the universal form (MDAnalysis
  `_fit_to`: `X' = R(X - Xbar) + Xbar_ref`; Eigen `umeyama`). **Do not "fix" the Kabsch.**
- The bug is isolated to `TransformedConformation::smoothTransformSequence`
  (`src/model/TransformedConformation.cpp:282-325`): it averages the rotation quaternion AND
  the translation **independently** (`out.T = tSum / count`, line 321). A smoothed R paired
  with an independently-smoothed T is no longer a matched centroid-pinned pair, so at the 93 Å
  lever arm the centroid swings. (Measured: centroid `|c|` 93 Å → 106 Å, Δz ≈ 36 Å by frame 750.)
- Second, quieter instability: the reference is frame 0 (`referenceFrame_`), an arbitrary
  snapshot, so even unsmoothed fits jitter.

## The fix — two parts, both grounded in the prior art
### A. Centroid-pin the smoothed transform (kills the walk)
Smooth ONLY the rotation. Re-derive the translation from the smoothed rotation and the
per-frame centroid so the fit-set centroid lands on a CONSTANT reference point every frame:
`T'_f = cr - R'_f * cc_f`, where `cc_f` = current frame's fit-set centroid, `cr` = the
(constant) reference-structure fit-set centroid. NEVER average T in world space.
- This requires `smoothTransformSequence` to know `cc_f` (per frame) and `cr` (constant).
  Today they live only inside the Kabsch call. Restructure so the raw-transform stage stores,
  per frame, the matched `cc_f` (and `cr` once) alongside `R_f` — e.g. extend the per-frame
  cache record, or carry parallel vectors. Keep `atomPosition` as `R*raw + T` (unchanged); only
  how `T` is produced changes.
- INVARIANT: for ANY smoothing window (including 0), the fit-set centroid maps exactly to `cr`.
  A quick self-check the harness can assert: centroid of the fit-set transformed positions is
  frame-invariant to < 1e-6 Å.

### B. Iterative converged-mean reference (kills the jitter) — MDAnalysis `iterative_average`
Replace the frame-0 reference with the iterative mean structure over the whole trajectory,
computed on the fit set (the same atoms the mode fits on: backbone subset for `FitSubset`, all
atoms for `FitReference`). Algorithm (cite MDAnalysis `iterative_average`,
https://docs.mdanalysis.org/stable/documentation_pages/analysis/align.html lines 513-529):
```
ref = fitset coords at frame 0, demeaned (centroid subtracted)
for iter in 0..maxIter:                      # maxIter ~ 20 (MDAnalysis default 100)
    accum = zeros(nFit, 3)
    for f in all frames:
        cur = fitset coords at f, demeaned
        R   = Kabsch_rotation(cur -> ref)    # rotation only (both already centered)
        accum += R * cur
    avg = accum / nFrames
    d   = rms(avg - ref)                      # RMSD between successive references
    ref = avg
    if d < eps: break                         # eps ~ 1e-4 Å (MDAnalysis 1e-6)
referencePositions_ = ref                      # the converged mean, in the fit-set order
cr = a CONSTANT world anchor for the centroid (frame-0 fit-set centroid is fine; document it)
```
- The converged mean replaces what `setMode` currently caches as `referencePositions_`
  (`src/model/TransformedConformation.cpp:143-155`). Per-frame `computeRawTransform`
  (`:205-254`) then Kabsch-fits each frame's fit set to this mean (unchanged in shape).
- Cost: ~20 iters × 751 frames × cheap Kabsch. Precompute at `setMode` (synchronous is OK if
  it's tens of ms; if it's clearly heavier, log the timing — do NOT silently stall the GUI).
- Both modes use the same machinery; the only difference stays the fit set.

### Smoothing default + knob
- Change `kDefaultStabilisationWindow` from 8 to **0** (`TransformedConformation.h:60`). No mature
  tool smooths the alignment transform; the stable mean reference is the stabiliser. KEEP the
  `setStabilisationWindow` setter + `POST /transform/smoothing` route as an optional cosmetic
  knob — but it now smooths only R and re-derives T (part A), so a non-zero window can never walk
  the molecule again.

## Truth bar (must hold)
Every per-frame transform is ONE rigid motion (a single R, plus the centroid map cc_f → cr)
applied to ALL atoms uniformly. That preserves all intra-molecular distances/angles by
construction, so real internal/domain motion survives — only global rigid-body motion is
removed. Do not introduce any per-atom or per-region transform.

## Out of scope for v1 (note, do not build)
Robust/weighted fitting for the all-atom "with give" gyroscope (Theseus-style ML downweighting
of variable regions; PyMOL/ChimeraX iterative outlier rejection). It's the documented cure for
all-atom precession and is a likely v2, but v1 ships iterative-mean + centroid-pin + stable-core
(backbone) only. Leave a one-line `// v2:` marker where weighting would slot in.

## Output
- Build green: `cmake --build --preset linux-rwdi --target h5reader -j$(nproc)` (or
  `cmake --build .../build/linux-rwdi --target h5reader`). Do NOT launch. Do NOT run git.
- Append a `notes/BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md` line: what landed (iterative-mean
  reference, centroid-pinned smoothing, default window 0), the convergence params used, and any
  REST/UI changes.
- Note anything in the code that contradicts this brief. If `cc_f`/`cr` plumbing forces a
  structural change you're unsure about, STOP and write the question rather than guess.
