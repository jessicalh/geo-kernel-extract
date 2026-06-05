# Codex — temporal stabilisation of the alignment transform (both stabilisation modes)

## (Standard qt6-cpp skill + VTK-docs-license preamble is prepended above — honour it.)

## Why
The two display-stabilisation modes (locked backbone = `FitSubset(BackboneSubset)`, with-give =
`FitReference` all-atom) each compute their alignment transform **per frame, independently** — `atomPosition()`
Kabsch-fits the current frame's subset to a cached reference and returns `R_f·raw + T_f`
(`src/model/TransformedConformation.cpp:51`, fit at `:157`). There is **no temporal coherence**, so the
per-frame rotation jitters/precesses over playback — the "Kabsch gyroscope." The fit already removes real
tumbling; what's left is alignment *jitter*. The lead's fix (correct): **smooth the transform sequence over
time**, then apply the smoothed transform to the raw atoms — orientation moves smoothly, every atom's raw
internal motion is preserved ("the give"). **Both modes get it; the only difference between them stays the fit
target.**

## Read first
- `TransformedConformation` (`src/model/TransformedConformation.{h,cpp}`): `Mode` = `FitReference` / `FitSubset`
  (`:58`); `atomPosition(frame, atom)` computes/caches the per-frame `(R,T)` and returns `R·raw+T` (`:51`); the
  Kabsch is `KabschFit`/`ComputeSubsetTransform` (`src/app/FitTargetMath.h`); reference positions cached from
  frame 0. Eigen is available (the math uses it).
- `frameCount()`, `atomPosition` raw source (the inner conformation's positions), and where `transformChanged`
  is emitted + consumed (`ReaderMainWindow.cpp` refresh).

## Deliverable — temporally-smoothed alignment
1. **Per-frame transforms across the trajectory.** For the active mode + subset, compute the per-frame Kabsch
   `(R_f, T_f)` mapping each frame's subset onto the reference-frame-0 subset, for all `f in [0, frameCount)`.
   (The trajectory is fully loaded; this is N cheap fits. Do it on mode/reference change; cache it.)
2. **Smooth the sequence over a symmetric window `W`** (tunable; default e.g. `W=8` frames → 17-frame window):
   - **Rotation** via `Eigen::Quaterniond`: `q_f` from `R_f`; for each centre `f`, average
     `q_k` over `k in [f-W, f+W]` (clamp at ends) with **double-cover sign alignment** (flip `q_k` if
     `q_k.dot(q_ref) < 0`, ref = the centre's quaternion), then **normalise** → `q'_f` → `R'_f`. (A windowed
     normalised quaternion mean is fine for small windows; if you prefer, an incremental SLERP smoothing is
     acceptable — your call, but it MUST handle the sign double-cover and renormalise.)
   - **Translation**: `T'_f = mean(T_k)` over the same clamped window.
3. **Apply the smoothed transform**: `atomPosition(f, atom) = R'_f · raw_f + T'_f`. Raw internal motion is
   untouched (still the raw frame's relative positions); only the rigid alignment is smoothed.
4. **The window is a tunable knob** so we can dial it visually:
   - a setter on `TransformedConformation` (e.g. `setStabilisationWindow(int halfWidth)`), recompute the
     smoothed sequence + emit `transformChanged`;
   - a REST route **`POST /transform/smoothing {"window": int}`** (and report it in `GET /transform`) so it can
     be tuned headlessly. `window=0` = no smoothing (current per-frame behaviour, for A/B).
5. Keep it correct for **both modes** (re-smooth when the mode or reference changes) and for the degenerate
   guards already in `KabschFit` (rank/determinant; a degenerate frame's transform should not poison the
   window — fall back to identity-rotation/centroid for that frame as today).

## Watch
- Quaternion double-cover (sign) is the classic bug — get it right or the smoothing flips/jerks.
- Don't smooth the *atom positions* (that would smear real motion) — smooth only the alignment `(R,T)`.
- Don't block the GUI thread on a huge recompute mid-playback; precompute on mode/window/reference change.
- The camera Kabsch path (Subset mode) shares `FitTargetMath` but is being CUT later — do not touch it here.

## Output
Build green (`cmake --build .../build/linux-rwdi --target h5reader -j$(nproc)`). Do NOT launch, do NOT git.
Append a `notes/BASIC_SCREEN_STATE_AND_PLAN_2026-06-05.md` line (what landed + the window default + the REST
route). Note any code that contradicts this brief.
